//
// Created by grenier on 28/05/24.
//

#ifndef ASTEX_CCVT_APPLICATION_H
#define ASTEX_CCVT_APPLICATION_H

#include "ccvt.h"

// Forward declare
struct GLFWwindow;

struct point_info{
    float _x;
    float _y;
    float _w;
    float _sel;

    point_info(float x, float y, float w):_x(x), _y(y), _w(w), _sel(0.){};
};

struct selection{
    bool active = false;
    int id = -1;
};

class ccvt_application {
public:
    // A function called only once at the beginning. Returns false is init failed.
    bool onInit();

    // A function called only once at the very end.
    void onFinish();

    // A function that tells if the application is still running.
    bool isRunning();

    // A function called at each frame, guaranteed never to be called before `onInit`.
    void onFrame();

private:

    bool initGui();
    void terminateGui();
    void updateGui();

    void displayPoints();
    void updatePoints();

    void onMouseButton(int button, int action, int mods);
    void onMouseMove(double xpos, double ypos);
    void addPoint(float xPos, float yPos);

private:
    GLFWwindow* m_window_H = nullptr;
    int m_width_H = 400;
    int m_height_H = 400;

    GLFWwindow* m_window_N1 = nullptr;
    GLFWwindow* m_window_N2 = nullptr;
    int m_width_N = 600;
    int m_height_N = 600;

    GLFWwindow* m_window_T = nullptr;
    int m_width_T = 800;
    int m_height_T = 800;

    unsigned int m_PointsShaderProgram;
    unsigned int m_PointsVAO;
    unsigned int m_PointsVBO;

    unsigned int m_ColorShaderProgram;
    unsigned int m_ColorVAO;

    unsigned int m_NoiseShaderProgram;
    unsigned int m_NoiseVAO;

    unsigned int m_CompositionShaderProgram;
    unsigned int m_CompositionVAO;

    unsigned int m_fbo_N1;
    unsigned int m_texture_N1;

    unsigned int m_fbo_N2;
    unsigned int m_texture_N2;


    CCVT m_ccvt;
    static const unsigned int m_MaxPointsNb = 12;
//    std::vector<Point> m_ccvtPoints{Point(0.81, 0.65),
//                                      Point(0.15, 0.25),
//                                      Point(0.52, 0.91),
//                                      Point(0.65, 0.28),
//                                      Point(0.12, 0.72),
//                                      Point(0.45, 0.48)};

//    float m_points[4*m_MaxPointsNb] = {0.81, 0.65, -0.0344047, 0.,
//                                       0.15, 0.25, 0.0645013, 0.,
//                                       0.52, 0.91 ,0.019120, 0.,
//                                       0.65, 0.28, -0.036392, 0.,
//                                       0.12, 0.72, 0.06444, 0.,
//                                       0.45, 0.48, -0.0390251, 0.};

//    std::vector<float> m_points = {0.81, 0.65, -0.0344047, 0.,
//                                       0.15, 0.25, 0.0645013, 0.,
//                                       0.52, 0.91 ,0.019120, 0.,
//                                       0.65, 0.28, -0.036392, 0.,
//                                       0.12, 0.72, 0.06444, 0.,
//                                       0.45, 0.48, -0.0390251, 0.};

    std::vector<point_info> m_points = {point_info{0.81, 0.65, 0.},
                                        point_info{0.15, 0.25, 0.},
                                        point_info{0.52, 0.91, 0.},
                                        point_info{0.65, 0.28, 0.},
                                        point_info{0.12, 0.72, 0.},
                                        point_info{0.45, 0.48, 0.}};

    selection m_selected;

    float m_F1Min = 10.;
    float m_F1Max = 10.;
    float m_Or1Min = 0.;
    float m_Or1Max = 0.;

    float m_F2Min = 12.;
    float m_F2Max = 12.;
    float m_Or2Min = 0.;
    float m_Or2Max = 3.;
};


#endif //ASTEX_CCVT_APPLICATION_H
