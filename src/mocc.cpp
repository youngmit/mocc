#include <iostream>
#include <iomanip>
#include <exception>
#include <sstream>

#include "global_config.hpp"
#include "files.hpp"
#include "input_proc.hpp"
#include "core_mesh.hpp"
#include "error.hpp"
#include "h5file.hpp"
#include "transport_sweeper.hpp"


using std::cout;
using std::cin;
using std::endl;

using namespace mocc;

/**
 * \mainpage MOCC
 *
 * Style
 * =====
 * For the most part, we are following the <a
 * href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.html">Google
 * C++ style guide</a>. Some highlights are below.
 *
 * Class Names
 * -----------
 * All class names should start with a capitol letter, and the first letter of
 * each subsequent word should be captiolized. Where one class is a simple or
 * subtle extension of another class or concept (e.g. MoCSweeper_2D3D), use an
 * underscore to separate the broad concept from the specialization. Most
 * extensions of classes should be written as normal (e.g. MoCSweeper as an
 * extension of TransportSweeper).
 *
 * File Naming
 * -----------
 * Most classes should have their own file, with the exception of small utility
 * classes which may make sense to be packaged together with a larger, more
 * standalone class. Class filenames should share the same name of the class
 * that they contain, be all lowercase, with words separated by underscores.
 *
 * Conventions
 * ===========
 * "A place for everything, and everything in its place."
 *
 * MOCC makes use of several conventions throughout the code in order to make
 * certain difficult aspects of medium- to large-scale code projects obvious and
 * simple.
 *
 * Constness
 * ---------
 * Be sure to understand the concept of constness. A few cood references can be
 * found <a href="http://www.cprogramming.com/tutorial/const_correctness.html">
 * here</a> and <a
 * href="http://duramecho.com/ComputerInformation/WhyHowCppConst.html">here</a>.
 * In general, make as many references, pointers, and class methods const as
 * possible, as it clearly specifies intent, protects the sanity and
 * ease-of-understanding of the data flow, and permits more optimization for the
 * compiler.
 *
 * Pointers and Ownership
 * ----------------------
 * C/C++ (especially C) does not provide a formal, language-level concept of
 * memory ownership; its kind of the Wild West when it comes to who should
 * destroy and deallocate what.
 *
 * The concept of "ownership" is simple in concept, but potentially a horrendous
 * mess if not approached with discipline. Ownership essentially means the scope
 * or program unit that is responsible for the appropriate maintenance of a
 * chunk of memory. For a scope or object to "own" a piece of memory means that
 * that scope is responsible for the proper destruction of that memory. If
 * everything is done properly, one should never encounter memory leaks or
 * double frees or the like. In practice, this is usually easier said than done,
 * as it is often not readily clear who the owner of an object (or rather a
 * reference to an object) is.
 *
 * To keep things reasonable, we will make heavy and deliberate use of the newer
 * memory management features of the C++ standard library in the \c memory
 * header, namely \c std::unique_ptr<T> and \c std::shared_ptr<T>. Raw pointers
 * are still used throughout, but their use is restricted to specific
 * circumstances. Here I will describe where each is used and why:
 *
 * ### <tt>shared_ptr<T></tt> ######
 * The shared pointer is used only when there is actual need for shared
 * ownership of memory (which in the case of MOCC should be rare). A shared
 * pointer contains an internal count of the number of places that are pointing
 * to the location that is wrapped by the shared pointer.  When any instance of
 * the shared pointer goes out of scope, the count is decremented, and whenever
 * a new shared pointer is associated with the location, the count is
 * incremented. The last shared pointer to the location pointed to is
 * responsible for destroying/deallocating it.
 *
 * This approach, while very flexible, is the least expressive in terms of what
 * it tells the developer about the ownership of the object being pointed to.
 * Indeed, it practically tells the developer nothing at all; potentially any
 * location that maintains an instance of the shared pointer could the the
 * ultimate custodian of the memory being pointed to, making guarantees about
 * object lifetime hard to determine.
 *
 * Due to the ambiguity built into \c shared_pointer<T>, it should only be used
 * in situations where the lifetime of the object is truely ambiguous and
 * difficult to determine at compile time. 
 *
 * \todo provide an example.
 *
 * ### <tt>unique_ptr<T></tt> ######
 * The unique pointer is much simpler than the shared pointer, in that it does
 * not maintain a count of references. It assumes that it is the sole owner
 * of the memory being pointed to (hence the name). The primary functionality
 * provided by the unique pointer is that its destructor <tt>delete</tt>s the
 * location with which it is associated, potentially calling that object's
 * destructor in the process.
 *
 * So the unique pointer provides a couple of things:
 *  - A guarantee to the developer that the lifetime of the object being pointed
 *  to is that of the unique pointer itself (and commonly, by extension, the
 *  lifetime of an object of which the unique pointer is an attribute).
 *  - The convenience of not having worry about explicitly <tt>delete</tt>ing
 *  the pointer, since it happens automatically when the unique pointer falls
 *  out of scope. Humans shouldn't have to remember such things; its below us.
 *
 * ### Raw pointers (<tt>T*</tt>) ######
 * Raw pointers are used to provide a reference to an object or data in which no
 * ownership is held by the pointer itself. Raw pointers are used in places
 * where a reference to an object is needed temporarily, and it is assumed that
 * ownership lies elsewhere. 
 *
 * \c new and \c delete should almost never be called on a raw pointer.
 *
 * Raw pointers are used sometimes as arguments to functions, which operate on
 * or otherwise use a reference to an object temporarily. They are also used
 * sometimes as data members of classes, instances of which are known not to
 * outlive the shared or unique pointer with which the raw pointer is
 * associated.
 *
 * An example of this are in the Pin types, which contain raw pointers to their
 * associated PinMesh.  Since the mocc::Pin and mocc::PinMesh objects are stored
 * in collections on the mocc::CoreMesh, we pretty much know that they have the
 * same lifetime, so using a raw pointer in the Pin class makes it clear that
 * the Pin does not own, and is not responsible for destroying the PinMesh
 * associated with it. CoreMesh contains a collection of unique pointers to
 * PinMesh objects, and is therefore the proper owner of them. The only
 * potentially scary thing here is that we need to know at the Pin level that
 * our lifetime will not extend past that of the PinMesh we are pointing to,
 * which is easy to enforce. Were that not the case, it might make more sense to
 * use shared pointers all around.
 *
 *
 * Coordinate Systems
 * ------------------
 * When working with the mocc::CoreMesh class, it is important to have a good
 * understanding of the various coordinate systems used. To facilitate a
 * simple-to-understand and highly-structured ray tracing approach, multiple
 * coordinate systems are used to conform to the heirarchical nature by which
 * geometry is specified in the first place (Pin -> Lattice -> Assembly ->
 * Core). While this makes ray tracing less of a mess in general, withouht an
 * understanding of the coordinate systems, it is likely harder to understand.
 *
 * Here are the different coordinate systems that will be encountered:
 *
 * ### Global (Core) Coordinates ####
 * Global coordinates are measured from the bottom, southwest corner of the
 * global spatial domain.
 *
 * ### Lattice Coordinates ####
 * Lattice coordinates are technically 2-D, and are measured from the southwest
 * corner of the lattice.
 *
 * ### Pin Coordinates ####
 * Pin coordinates are also 2-D, and mesasured from the dead center of the pin
 * geometry. This is so that the most common case, a cylindrical pin, centered
 * in the pin cell is easiest to trace. When encountering code that provides the
 * location of a pin origin from in core/global coordinates, this specifies the
 * origin of the pin coordinate system, translated to its location in real
 * space, then converted to the global coordinate scheme. This is done to make
 * conversion of global-coordinate points to pin-local points easy (refer to
 * mocc::RayData::RayData() and mocc::CoreMesh::get_pinmesh() to see how this is
 * done.)
 */

// The global top-level solver
SP_Solver_t solver;

// Global core mesh
SP_CoreMesh_t mesh;

int main(int argc, char* argv[]){
	// Make sure we have an input file
	if(argc < 2){
		Error("No input file specified!");
	}
	
	// Spin up the log file. For now, just use the name of the input file.
	StartLogFile(argv[1]);
	
	LogFile << "Welcome to " << PROG_NAME << "!" << std::endl << std::endl;

	// Parse the input file
	InputProc inProc(argv[1]);

    // Get an SP to the core mesh
    mesh = inProc.core_mesh();

    // Pull a shared pointer to the top-level solver and make it go
    solver = inProc.solver();
    solver->solve();

    // Output stuff
    H5File outfile("out.h5");
    solver->output( outfile );

    StopLogFile();
}
