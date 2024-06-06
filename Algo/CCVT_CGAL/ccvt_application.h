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

struct cell_info{
    float _cap;
    float _r = 0.;
    float _g = 0.;
    float _b = 0.;

    cell_info(float cap, float r, float g, float b):_cap(cap), _r(r), _g(g), _b(b){};
};

struct selection{
    bool active = false;
    bool toogle = false;
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

    void attacheCCVT(CCVT &ccvt){
        m_ccvt = ccvt;
    }

private:
    bool initGui(GLFWwindow* &window);
    void terminateGui();
    void updateGui();

    void displayPoints();
    void updatePoints();

    void onMouseButton(int button, int action, int mods);
    void onMouseMove(double xpos, double ypos);
    void addPoint(float xPos, float yPos);
    void deletePoint(int id);

    void getCCVTcells();
    void optimizeCCVT();
    void updateCCVT();

private:
    GLFWwindow* m_window_H = nullptr;
    int m_width_H = 400;
    int m_height_H = 400;

    GLFWwindow* m_window_N1 = nullptr;
    GLFWwindow* m_window_N2 = nullptr;
    int m_width_N = 400;
    int m_height_N = 400;

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

    unsigned int m_fbo_H;
    unsigned int m_texture_H;


    CCVT m_ccvt;
    static const unsigned int m_MaxPointsNb = 30;
    unsigned int m_CurrentPointsNb;

    std::vector<point_info> m_points;
    std::vector<cell_info> m_cells;

    selection m_selected;


    float m_F1Princ = 10.;
    float m_F1Spread = 0.;
    float m_Or1Princ = 2.;
    float m_Or1Spread = 0.;

    float m_F2Princ = 16.;
    float m_F2Spread = 0.;
    float m_Or2Princ = 1.;
    float m_Or2Spread = 3.;
};


#endif //ASTEX_CCVT_APPLICATION_H
