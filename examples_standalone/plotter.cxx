#include "GamPolaPlotter/GGampolaPlotter.h"

using namespace Gamapola;

int main()
{
    auto&& plotter = std::make_shared<GGampolaPlotter>("../additionalFiles/plot.img");
    plotter->plot("neutral_100k.root", 0, "neutral_100k_plot.root");
    return 0;
}


