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

struct histo_info{
    int _count;
    float _obj;
    float _r;
    float _g;
    float _b;

    histo_info(float r, float g, float b): _count(0), _obj(0.), _r(r), _g(g), _b(b){};
    void incr(){_count++;}
};

struct selection{
    bool active = false;
//    bool toogle = false;
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
    void drawNoise(unsigned int &fbo, float width, float height, float freq_princ, float freq_spread, float omega_princ, float omega_spread, float seed);
    void drawCM();
    void drawComposition();

    bool initGui(GLFWwindow* &window);
    void terminateGui();
    void updateGui();
    void cellPopup(int id);
//
//    void displayPoints();
//    void updatePoints();
//
//    void onMouseButton(int button, int action, int mods);
//    void onMouseMove(double xpos, double ypos);
    void addPoint(float xPos, float yPos, float r, float g, float b);
    void insertPoint(int vecPos, float xPos, float yPos, float r, float g, float b);
    void deletePoint(int id);
    void normilizeCap();
    void equalizeCap();

    void getCCVTcells();
    void optimizeCCVT();
    void updateCCVT();

    void saveTexture(unsigned  int fbo_id, int width, int height, const std::string &filename);
//    void savePPM(const std::string& filename, unsigned int* data, int width, int height);


    void computeStatistiques(unsigned int fbo_id, int width, int height, float &mean, float &var);
    float computeMean(std::vector<float> data);
    float computeSquareMean(std::vector<float> data);
    void computeProportions();



private:
//    GLFWwindow* m_window_H = nullptr;
    int m_width_H = 400;
    int m_height_H = 400;

//    GLFWwindow* m_window_N1 = nullptr;
//    GLFWwindow* m_window_N2 = nullptr;
    int m_width_N = 400;
    int m_height_N = 400;

    GLFWwindow* m_window_T = nullptr;
    int m_width_T = 1600;
    int m_height_T = 900; // 900

//    unsigned int m_PointsShaderProgram;
//    unsigned int m_PointsVAO;
//    unsigned int m_PointsVBO;

    unsigned int m_ColorShaderProgram;
    unsigned int m_ColorVAO;

    unsigned int m_NoiseShaderProgram;
    unsigned int m_NoiseVAO;

    unsigned int m_CompositionShaderProgram;
    unsigned int m_CompositionVAO;

    unsigned int m_fbo_N1;
    unsigned int m_texture_N1;
    unsigned int m_fbo_N1_ui;
    unsigned int m_texture_N1_ui;

    unsigned int m_fbo_N2;
    unsigned int m_texture_N2;
    unsigned int m_fbo_N2_ui;
    unsigned int m_texture_N2_ui;

    unsigned int m_fbo_H;
    unsigned int m_texture_H;
    unsigned int m_fbo_T;
    unsigned int m_texture_T;

    bool m_n1_changed = false;
    bool m_n2_changed = false;

    float m_mean_N1;
    float m_mean_N2;
    float m_var_N1;
    float m_var_N2;



    CCVT m_ccvt;
    static const unsigned int m_MaxPointsNb = 30; // change in composition shader and color shader too
    unsigned int m_CurrentPointsNb;

    std::vector<point_info> m_points;
    std::vector<cell_info> m_cells;
    std::vector<histo_info> m_histo;

    selection m_selected;

    std::string m_infoBuffer = "";


    float m_F1Princ = 16.;
    float m_F1Spread = 0.;
    float m_Or1Princ = 2.;
    float m_Or1Spread = 0.;

    float m_F2Princ = 20.;
    float m_F2Spread = 0.;
    float m_Or2Princ = 1.;
    float m_Or2Spread = 3.;
};


#endif //ASTEX_CCVT_APPLICATION_H
