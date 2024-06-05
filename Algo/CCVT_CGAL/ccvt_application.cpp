//
// Created by grenier on 28/05/24.
//

#include "ccvt_application.h"

#include <GL/glew.h>
#include <GL/gl.h>
#include <GLFW/glfw3.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <iostream>
#include <fstream>
#include <sstream>

void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_compilation(unsigned int shader){
    // checking the shader compilation
    int success;
    char infoLog[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if(!success){
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        std::cerr<<"shader compilation failed\n"<<infoLog<<std::endl;
    }
}

void check_link(unsigned int shader){
    // checking the linking of program
    int success;
    char infoLog[512];
    glGetProgramiv(shader, GL_LINK_STATUS, &success);
    if(!success){
        glGetProgramInfoLog(shader, 512, nullptr, infoLog);
        std::cerr<<"shader linking failed\n"<<infoLog<<std::endl;
    }
}

std::string load_shader(const std::string& file_path){
    std::cout<<"loading shader "<<file_path<<"..."<<std::endl;
    std::string data;
    std::ifstream dataStream(file_path);

    if(dataStream.is_open()){
        std::stringstream sstr;
        sstr << dataStream.rdbuf();
        data = sstr.str();
    }else{
        std::cerr<<"can't open "<<file_path<<" check the working directory"<<std::endl;
    }
    return data;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool window_creation(GLFWwindow* &window, const char* name, int width, int height){
    std::cout<<"creation window "<<name<<"..."<<std::endl;

    // fenêtre carte H
    window = glfwCreateWindow(width, height, name, NULL, NULL);
    if (!window) {
        std::cerr << "Could not open window!" << std::endl;
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window);


    // initialising GLEW
    glewExperimental = GL_TRUE;
    auto err = glewInit();
    if(err != GLEW_OK){
        std::cerr<<"Fail to initialize GLEW"<<std::endl;
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
        return false;
    }

    return true;
}


void display_quad(unsigned int & VAO){
    // vertex
    float vertices[] = {1.f, 1.f,
                        1.f, -1.f,
                        -1.f, -1.f,
                        -1.f, 1.f};
    unsigned int indices[] = {0, 1, 3,
                              1, 2, 3};

    // Vertex Array object
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO); // 1. bind VAO

    // Vertex buffer object
    unsigned int VBO;
    glGenBuffers(1, &VBO); // id du buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO); // 2. copy the vertices in a buffer
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // Elements buffer object
    unsigned int EBO;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);// 3. copy the index in a buffer
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);


    glVertexAttribPointer( // 4. set the vertex attributes pointer
            0, // location
            2, // size of attribute (vec2)
            GL_FLOAT, // data type
            GL_FALSE, // data need to bee normalized ?
            2*sizeof(float), // distance between consecutive vertex
            (void*)0); // offset
    glEnableVertexAttribArray(0);
}



void shader_program(unsigned int & shaderProgram, std::string frag_shader, std::string vert_shader){
    // vertex shader
    std::string vertexShaderSourceString = load_shader(vert_shader);
    const char* vertexShaderSource = vertexShaderSourceString.c_str();
    unsigned int vertexShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);
    check_compilation(vertexShader);

    // fragment shader
    std::string fragmentShaderSourceString = load_shader(frag_shader);
    const char* fragmentShaderSource = fragmentShaderSourceString.c_str();
    unsigned int fragmentShader;
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);
    check_compilation(fragmentShader);

    // shader program (linkage)
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    check_link(shaderProgram);

    // using shaders program (and delete shaders)
    glUseProgram(shaderProgram);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
}


bool setFBO(unsigned int &fbo, unsigned int &texture, int width, int height){
    // fbo noise 1
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE){
        std::cerr<<"erreur frame buffer"<< std::endl;
        return false;
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return true;
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ccvt_application::onInit() {
    std::cout<<"Initialisation..."<<std::endl;
    glfwSetErrorCallback(error_callback);

    // information initiale du ccvt
    getCCVTcells();


    // initialisation de glfw
    auto err_glfw = glfwInit();
    if (err_glfw != GLFW_TRUE) {
        std::cerr << "Could not initialize GLFW!" << std::endl;

        return false;
    }


    // création de fenêtres
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);


    // fenêtre carte H
    window_creation(m_window_H, "CCVT application", m_width_H, m_height_H);

    // Set the user pointer to be "this"
    glfwSetWindowUserPointer(m_window_H, this);
    // mouse button
    glfwSetMouseButtonCallback(m_window_H, [](GLFWwindow* window, int button, int action, int mods) {
        auto that = reinterpret_cast<ccvt_application*>(glfwGetWindowUserPointer(window));
        if (that != nullptr) that->onMouseButton(button, action, mods);
    });
    // mouse move
    glfwSetCursorPosCallback(m_window_H, [](GLFWwindow* window, double xpos, double ypos) {
        auto that = reinterpret_cast<ccvt_application*>(glfwGetWindowUserPointer(window));
        if (that != nullptr) that->onMouseMove(xpos, ypos);
    });

    shader_program(m_PointsShaderProgram, "shaders/points_shader.frag", "shaders/points_shader.vert");
    displayPoints();
    shader_program(m_ColorShaderProgram, "shaders/color_shader.frag", "shaders/identity_shader.vert");
    display_quad(m_ColorVAO);







    // fenêtre résultat
    window_creation(m_window_T, "composition", m_width_T, m_height_T);

    shader_program(m_NoiseShaderProgram, "shaders/noise_shader.frag", "shaders/identity_shader.vert");
    display_quad(m_NoiseVAO);

    shader_program(m_CompositionShaderProgram, "shaders/composition_shader.frag", "shaders/identity_shader.vert");
    display_quad(m_CompositionVAO);


    // fbo noise 1
    if(!setFBO(m_fbo_N1, m_texture_N1, m_width_T, m_height_T)){
        return false;
    }

    // fbo noise 2
    if(!setFBO(m_fbo_N2, m_texture_N2, m_width_T, m_height_T)){
        return false;
    }

    initGui(m_window_T);



    return true;
}



void ccvt_application::onFinish() {
    terminateGui();

    glDeleteVertexArrays(1, &m_PointsVAO);
    glDeleteProgram(m_PointsShaderProgram);

    glDeleteVertexArrays(1, &m_ColorVAO);
    glDeleteProgram(m_ColorShaderProgram);
//
    glDeleteVertexArrays(1, &m_NoiseVAO);
    glDeleteProgram(m_NoiseShaderProgram);

    glDeleteVertexArrays(1, &m_CompositionVAO);
    glDeleteProgram(m_CompositionShaderProgram);

    glfwDestroyWindow(m_window_H);
    glfwDestroyWindow(m_window_T);
    glfwTerminate();
}



bool ccvt_application::isRunning() {
    // close the window when "escape" is press
    bool is_escape_pressed = glfwGetKey(m_window_H, GLFW_KEY_ESCAPE) == GLFW_PRESS
                            or glfwGetKey(m_window_T, GLFW_KEY_ESCAPE) == GLFW_PRESS;

    if(is_escape_pressed){
        glfwSetWindowShouldClose(m_window_H, true);
        glfwSetWindowShouldClose(m_window_T, true);
    }

    return !glfwWindowShouldClose(m_window_H);
}





void ccvt_application::onFrame() {




    ////////////////////////////////////////////////////////////////////
    glfwMakeContextCurrent(m_window_H);
    glViewport(0, 0, m_width_H, m_height_H); // (lower_left_x, lower_left_y, width, height)

    // check and call events
    glfwPollEvents();

    // clear color
    glClearColor(0.4f, 0.6f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // colors
    glUseProgram(m_ColorShaderProgram);

    // unifoms
    unsigned int UBO_points; // uniform buffer
    glGenBuffers(1, &UBO_points);
    glBindBuffer(GL_UNIFORM_BUFFER, UBO_points);
    glBufferData(GL_UNIFORM_BUFFER, m_MaxPointsNb*sizeof(point_info), m_points.data(), GL_STREAM_DRAW);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, UBO_points, 0, m_MaxPointsNb * sizeof(point_info));

    unsigned int UBO_colors; // uniform buffer
    glGenBuffers(1, &UBO_colors);
    glBindBuffer(GL_UNIFORM_BUFFER, UBO_colors);
    glBufferData(GL_UNIFORM_BUFFER, m_MaxPointsNb*sizeof(cell_info), m_cells.data(), GL_STREAM_DRAW);
    glBindBufferRange(GL_UNIFORM_BUFFER, 1, UBO_colors, 0, m_MaxPointsNb * sizeof(cell_info));

    glUniformBlockBinding(m_ColorShaderProgram, glGetUniformBlockIndex(m_ColorShaderProgram, "uPoints"), 0);
    glUniformBlockBinding(m_ColorShaderProgram, glGetUniformBlockIndex(m_ColorShaderProgram, "uColors"), 1);
    glUniform2f(glGetUniformLocation(m_ColorShaderProgram, "uRes"), m_width_H, m_height_H);
    glUniform1i(glGetUniformLocation(m_ColorShaderProgram, "uPointsNb"), m_CurrentPointsNb);


    glBindVertexArray(m_ColorVAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glBindBuffer(GL_UNIFORM_BUFFER, 0);



    // points
    glUseProgram(m_PointsShaderProgram);
    glBindVertexArray(m_PointsVAO);
    glDrawArrays(GL_POINTS, // mode
                   0, // first
                   m_CurrentPointsNb); // count


    // swap the buffers
    glfwSwapBuffers(m_window_H);










    ////////////////////////////////////////////////////////////////////
    glfwMakeContextCurrent(m_window_T);
    glViewport(0, 0, m_width_T, m_height_T); // (lower_left_x, lower_left_y, width, height)


    // fbo noise 1
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_N1);

    // clear color
    glClearColor(0.6f, 0.6f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // we're not using the stencil buffer now
    glEnable(GL_DEPTH_TEST);

    glUseProgram(m_NoiseShaderProgram);

    // unifoms
    glUniform2f(glGetUniformLocation(m_NoiseShaderProgram, "uRes"), m_width_T, m_height_T);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uFprinc"), m_F1Princ);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uFspread"), m_F1Spread);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uOprinc"), m_Or1Princ);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uOspread"), m_Or1Spread);

    glBindVertexArray(m_NoiseVAO);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);





    // fbo noise 2
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_N2);

    // clear color
    glClearColor(0.6f, 0.6f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // we're not using the stencil buffer now
    glEnable(GL_DEPTH_TEST);

    glUseProgram(m_NoiseShaderProgram);

    // unifoms
    glUniform2f(glGetUniformLocation(m_NoiseShaderProgram, "uRes"), m_width_T, m_height_T);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uFprinc"), m_F2Princ);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uFspread"), m_F2Spread);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uOprinc"), m_Or2Princ);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uOspread"), m_Or2Spread);

    glBindVertexArray(m_NoiseVAO);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);




    // composition
    glClearColor(0.4f, 0.6f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    glUseProgram(m_CompositionShaderProgram);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, m_texture_N1);
    glUniform1i(glGetUniformLocation(m_CompositionShaderProgram, "uTex_1"), 0);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, m_texture_N2);
    glUniform1i(glGetUniformLocation(m_CompositionShaderProgram, "uTex_2"), 1);


    // uniform buffer
    unsigned int UBO_H_points;
    glGenBuffers(1, &UBO_H_points);
    glBindBuffer(GL_UNIFORM_BUFFER, UBO_H_points);
    glBufferData(GL_UNIFORM_BUFFER, m_MaxPointsNb*sizeof(point_info), m_points.data(), GL_STREAM_DRAW);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, UBO_H_points, 0, m_MaxPointsNb * sizeof(point_info));

    unsigned int UBO_H_colors; // uniform buffer
    glGenBuffers(1, &UBO_H_colors);
    glBindBuffer(GL_UNIFORM_BUFFER, UBO_H_colors);
    glBufferData(GL_UNIFORM_BUFFER, m_MaxPointsNb*sizeof(cell_info), m_cells.data(), GL_STREAM_DRAW);
    glBindBufferRange(GL_UNIFORM_BUFFER, 1, UBO_H_colors, 0, m_MaxPointsNb * sizeof(cell_info));

    glUniformBlockBinding(m_CompositionShaderProgram, glGetUniformBlockIndex(m_CompositionShaderProgram, "uPoints"), 0);
    glUniformBlockBinding(m_CompositionShaderProgram, glGetUniformBlockIndex(m_CompositionShaderProgram, "uColors"), 1);
    glUniform2f(glGetUniformLocation(m_CompositionShaderProgram, "uRes"), m_width_T, m_height_T);
    glUniform1i(glGetUniformLocation(m_CompositionShaderProgram, "uPointsNb"), m_CurrentPointsNb);


    glBindVertexArray(m_CompositionVAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glBindBuffer(GL_UNIFORM_BUFFER, 0);


    updateGui();


    // swap the buffers
    glfwSwapBuffers(m_window_T);



}




////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ccvt_application::displayPoints() {
    glEnable(GL_PROGRAM_POINT_SIZE);


    // Vertex Array object
    glGenVertexArrays(1, &m_PointsVAO);

    // Vertex buffer object
    glGenBuffers(1, &m_PointsVBO);

    updatePoints();

}


void ccvt_application::updatePoints() {

//    for(auto pi:m_points){
//        std::cout<<pi._x<<", "<<pi._y<<std::endl;
//    }

    // Vertex Array object
    glBindVertexArray(m_PointsVAO); // 1. bind VAO

    // Vertex buffer object
    glBindBuffer(GL_ARRAY_BUFFER, m_PointsVBO); // 2. copy the vertices in a buffer
    glBufferData(GL_ARRAY_BUFFER, m_MaxPointsNb*sizeof(point_info), m_points.data(), GL_STREAM_DRAW);


    glVertexAttribPointer( // 4. set the vertex attributes pointer
            0, // location
            4, // size of attribute (vec4)
            GL_FLOAT, // data type
            GL_FALSE, // data need to bee normalized ?
            4*sizeof(float), // distance between consecutive vertex
            (void*)0); // offset
    glEnableVertexAttribArray(0);
}









////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ccvt_application::initGui(GLFWwindow* &window) {
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

// Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);          // Second param install_callback=true will install GLFW callbacks and chain to existing ones.
    ImGui_ImplOpenGL3_Init();

    return true;
}


void ccvt_application::terminateGui() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}


void ccvt_application::updateGui() {
    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    ImGui::ShowDemoWindow(); // Show demo window! :)

    ImGui::Begin("Hello, world!");
    ImGuiIO& io = ImGui::GetIO();
    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
    ImGui::End();


//    ImGui::Begin("color map");
//    bool optimize = ImGui::Button("volume optimisation");
//    if(optimize){
//        optimizeCCVT();
//    }
//    ImGui::End();



    ImGui::Begin("noise 1");
    ImGui::SliderFloat("F_0", &m_F1Princ, 4., 30.);
    ImGui::SliderFloat("F spread", &m_F1Spread, 0., 20.);
    ImGui::SliderAngle("O_0", &m_Or1Princ);
    ImGui::SliderAngle("O spread", &m_Or1Spread, 0., 180.);

    // we access the ImGui window size
    const float window_width_N1 = m_width_N;//ImGui::GetContentRegionAvail().x;
    const float window_height_N1 = m_height_N;//ImGui::GetContentRegionAvail().y;

    glViewport(0, 0, window_width_N1, window_height_N1);

    // we get the screen position of the window
    ImVec2 pos_N1 = ImGui::GetCursorScreenPos();

    // and here we can add our created texture as image to ImGui
    // unfortunately we need to use the cast to void* or I didn't find another way tbh
    ImGui::GetWindowDrawList()->AddImage(
            (void *)m_texture_N1,
            ImVec2(pos_N1.x, pos_N1.y),
            ImVec2(pos_N1.x + window_width_N1, pos_N1.y + window_height_N1),
            ImVec2(0, 1),
            ImVec2(1, 0)
    );
    ImGui::End();



    ImGui::Begin("noise 2");
    ImGui::SliderFloat("F_0", &m_F2Princ, 4., 30.);
    ImGui::SliderFloat("F spread", &m_F2Spread, 0., 20.);
    ImGui::SliderAngle("O_0", &m_Or2Princ);
    ImGui::SliderAngle("O spread", &m_Or2Spread, 0., 180.);

    // we access the ImGui window size
    const float window_width_N2 = m_width_N;//ImGui::GetContentRegionAvail().x;
    const float window_height_N2 = m_height_N;//ImGui::GetContentRegionAvail().y;

    glViewport(0, 0, window_width_N2, window_height_N2);

    // we get the screen position of the window
    ImVec2 pos_N2 = ImGui::GetCursorScreenPos();

    // and here we can add our created texture as image to ImGui
    // unfortunately we need to use the cast to void* or I didn't find another way tbh
    ImGui::GetWindowDrawList()->AddImage(
            (void *)m_texture_N2,
            ImVec2(pos_N2.x, pos_N2.y),
            ImVec2(pos_N2.x + window_width_N2, pos_N2.y + window_height_N2),
            ImVec2(0, 1),
            ImVec2(1, 0)
    );
    ImGui::End();


    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ccvt_application::onMouseButton(int button, int action, int mods) {
//    ImGuiIO& io = ImGui::GetIO();
//    if (io.WantCaptureMouse) {
//        // Don't rotate the camera if the mouse is already captured by an ImGui
//        // interaction at this frame.
//        return;
//    }


    if(button == GLFW_MOUSE_BUTTON_LEFT)
    {
        switch (action) {
            case GLFW_PRESS:
                m_selected.toogle = true;
                if(!m_selected.active){
                    double xpos;
                    double ypos;
                    glfwGetCursorPos(m_window_H, &xpos, &ypos);
                    addPoint(xpos / m_width_H, 1.- ypos / m_height_H);
                }
                break;
            case GLFW_RELEASE:
                m_selected.toogle = false;
                m_selected.active = false;
                m_selected.id = -1;

//                updateCCVT();
                break;
        }
    }

    if(button == GLFW_MOUSE_BUTTON_RIGHT and action == GLFW_PRESS)
    {
        optimizeCCVT();
    }
}

void ccvt_application::onMouseMove(double xpos, double ypos) {
    float xpos_normalized = xpos / m_width_H;
    float ypos_normalized = 1.- ypos / m_height_H;

//    std::cout<<m_selected.active<<" "<<m_selected.toogle<<" "<<m_selected.id<<std::endl;

    if(!m_selected.toogle){
        m_selected.active = false;
        m_selected.id = -1;
        for(auto p=m_points.begin(); p<m_points.end(); p++){
            float len = (xpos_normalized - (*p)._x)*(xpos_normalized - (*p)._x) + (ypos_normalized - (*p)._y)*(ypos_normalized - (*p)._y);

            if(len < .01){
                m_selected.active = true;
                m_selected.id = p-m_points.begin();
                (*p)._sel = 1.;
            }
            else{
                (*p)._sel = 0.;
            }
        }
    }

    if(m_selected.active and m_selected.toogle and m_selected.id != -1){
        m_points.at(m_selected.id)._x = xpos_normalized;
        m_points.at(m_selected.id)._y = ypos_normalized;
    }

    updatePoints();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ccvt_application::addPoint(float xPos, float yPos) {
//    std::cout<<"new point : "<<xPos<<", "<<yPos<<std::endl;

    float new_cap = 1./(m_cells.size()+1);

    float caps_sum = new_cap;
    for(auto ci:m_cells){
        caps_sum += ci._cap;
    }
    for(auto &ci:m_cells){
        ci._cap = ci._cap/caps_sum;
    }

    m_points.emplace_back(xPos, yPos, 0.);
    m_cells.emplace_back(new_cap/caps_sum, 1., 1., 1.);
    m_CurrentPointsNb = m_points.size();


    updatePoints();
}




void ccvt_application::getCCVTcells(){
    std::vector<FT> weights;
    std::vector<Point> points;
    std::vector<FT> proportions;
    std::vector<double> R;
    std::vector<double> G;
    std::vector<double> B;

    std::vector<point_info> tmp_points;
    std::vector<cell_info> tmp_cells;

    m_ccvt.collect_sites(points, weights);
    m_ccvt.get_proportion(proportions);
    m_ccvt.get_colors(R, G, B);

    for(int i=0; i<points.size(); i++){
        tmp_points.emplace_back(points.at(i).x()+0.5, points.at(i).y()+0.5, weights.at(i));
        tmp_cells.emplace_back(proportions.at(i), R.at(i), G.at(i), B.at(i));
    }

    m_points = tmp_points;
    m_cells = tmp_cells;
    m_CurrentPointsNb = m_points.size();
}




void ccvt_application::optimizeCCVT(){
    updateCCVT();
    std::cout<<"optimizing CCVT..."<<std::endl;
    FT stepX = 0.01; // pour les positions
    FT stepW = 0.1; // pour les poids
    FT epsilon = 1.;
    unsigned max_newton_iters = 500;
    unsigned max_iters = 500;

    int iter_opt = m_ccvt.optimize_H(stepW, stepX, max_newton_iters, epsilon, max_iters);

    std::cout<<iter_opt<<" itérations (max : "<<max_newton_iters<<")"<<std::endl;

    getCCVTcells();
    updatePoints();
}




void ccvt_application::updateCCVT(){
    std::cout<<"update CCVT..."<<std::endl;

    std::vector<Point> points;
    std::vector<double> weights;
    std::vector<double> caps;
    std::vector<double> R;
    std::vector<double> G;
    std::vector<double> B;

    for(auto pi:m_points){
        points.emplace_back(pi._x, pi._y);
        weights.emplace_back(pi._w);
    }

    for(auto ci:m_cells){
        caps.emplace_back(ci._cap);
        R.emplace_back(ci._r);
        G.emplace_back(ci._g);
        B.emplace_back(ci._b);
    }

    m_ccvt.set_custom_proportions(caps);
    m_ccvt.set_sites(points, weights);
    m_ccvt.set_colors(R, G, B);


//    std::cout<<"vérification"<<std::endl;
//    std::vector<FT> weights_t;
//    std::vector<Point> points_t;
//    std::vector<FT> proportions_t;
//    m_ccvt.collect_sites(points_t, weights_t);
//    m_ccvt.get_proportion(proportions_t);
//    for(int i=0; i<points_t.size(); i++){
//        std::cout<<points_t.at(i).x()+0.5<<", "<<points_t.at(i).y()+0.5<<", "<<weights_t.at(i)<<std::endl;
//    }

}