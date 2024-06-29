//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_TIMER_H
#define CCVT_TEST_CHARLINE_TIMER_H


// STL
#include <ctime>
#include <vector>
#include <string>
#include <iostream>

// local
#include "console_color.h"

class Timer
{
public:
    static
    void start_timer(std::vector<double>& timer,
                     const int color,
                     const std::string msg)
    {
        set_color(color);
        std::cout << msg << white << " ... ";
        timer.push_back(clock());
    }

    static
    double stop_timer(std::vector<double>& timer,
                    const int color,
                    const std::string msg = std::string())
    {
        double duration = time_duration(timer.back());
        timer.pop_back();
        set_color(color);
        std::cout << "done" << white << " (" << duration << " s) " << msg << std::endl;
        return duration;
    }

    static
    void set_color(const int color)
    {
        switch (color)
        {
            case COLOR_WHITE:
                std::cout << white;
                break;
            case COLOR_RED:
                std::cout << red;
                break;
            case COLOR_BLUE:
                std::cout << blue;
                break;
            case COLOR_GREEN:
                std::cout << green;
                break;
            case COLOR_YELLOW:
                std::cout << yellow;
                break;
        }
    }

    static
    double time_duration(const double init)
    {
        return (clock() - init) / CLOCKS_PER_SEC;
    }
};

#endif //CCVT_TEST_CHARLINE_TIMER_H
