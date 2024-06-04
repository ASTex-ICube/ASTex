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


    void getCCVTcells(){
        std::vector<FT> weights;
        std::vector<Point> points;
        std::vector<point_info> tmp_points;

        m_ccvt.collect_sites(points, weights);

        for(int i=0; i<points.size(); i++){
            tmp_points.push_back(point_info(points.at(i).x()+0.5, points.at(i).y()+0.5, weights.at(i)));
        }
        m_points = tmp_points;
    }

    void optimizeCCVT(){
        std::cout<<"optimizing CCVT..."<<std::endl;
        FT stepX = 0.01; // pour les positions
        FT stepW = 0.1; // pour les poids
        FT epsilon = 1.;
        unsigned max_newton_iters = 10;// 500;
        unsigned max_iters = 10;// 500;

        m_ccvt.optimize_H(stepW, stepX, max_newton_iters, epsilon, max_iters);
        getCCVTcells();
        updatePoints();
    }

    void updateCCVT(){
        std::cout<<"update CCVT..."<<std::endl;
        std::vector<Point> points;
        for(auto pi:m_points){
            points.push_back(Point{pi._x, pi._y});
        }
        m_ccvt.set_initial_sites(points);
    }

private:
    GLFWwindow* m_window_H = nullptr;
    int m_width_H = 600;
    int m_height_H = 600;

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

    unsigned int m_meanShaderProgram;

    unsigned int m_fbo_N1;
    unsigned int m_texture_N1;

    unsigned int m_fbo_N2;
    unsigned int m_texture_N2;


    CCVT m_ccvt;
    static const unsigned int m_MaxPointsNb = 12;
//    std::vector<point_info> m_points = {point_info{0.81, 0.65, 0.},
//                                        point_info{0.15, 0.25, 0.},
//                                        point_info{0.52, 0.91, 0.},
//                                        point_info{0.65, 0.28, 0.},
//                                        point_info{0.12, 0.72, 0.},
//                                        point_info{0.45, 0.48, 0.}};
    std::vector<point_info> m_points;

    selection m_selected;


    float m_F1Princ = 10.;
    float m_F1Spread = 0.;
    float m_Or1Princ = 1.;
    float m_Or1Spread = 0.;

    float m_F2Princ = 16.;
    float m_F2Spread = 0.;
    float m_Or2Princ = 1.;
    float m_Or2Spread = 3.;
};


#endif //ASTEX_CCVT_APPLICATION_H
