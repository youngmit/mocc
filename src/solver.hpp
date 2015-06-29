#pragma once

// This provides a virtual base type, which shall provide a solve() and step()
// method. At the highest level of the heirarchy, the driver calls solve() and
// that should invoke everything that is necessary.
class Solver{
    virtual void solve()=0;
};
