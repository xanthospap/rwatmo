#include "vmf3_grid_stream.hpp"

using namespace dso;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "Error. Usage: %s [VMF3GR GRID FILE]\n", argv[0]);
        return 1;
    }

    vmf3::GridVmf3Data grid;
    if (load_vmfgr3_grid_map(argv[1], &grid))
    {
        fprintf(stderr, "Error. Failed loading grid file %s\n", argv[1]);
        return 1;
    }

    return 0;
}