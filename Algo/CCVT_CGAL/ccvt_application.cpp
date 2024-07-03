//
// Created by grenier on 28/05/24.
//

#include "ccvt_application.h"
#include "timer.h"

#include <GL/glew.h>
#include <GL/gl.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "imgui/implot.h"

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
bool window_creation(GLFWwindow* &window, const char* name, int& width, int& height){
    std::cout<<"creation window "<<name<<"..."<<std::endl;

    // fenêtre carte H
    window = glfwCreateWindow(width, height, name, NULL, NULL);
    if (!window) {
        std::cerr << "Could not open window!" << std::endl;
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window);

    // int width_t, height_t;
    // glfwGetWindowSize(window, &width_t, &height_t);
    glfwGetFramebufferSize(window, &width, &height);
    // std::cout<<width_t<<", "<<height_t<<std::endl;
    // glfwSetWindowSize(window, width, height);

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





    // fenêtre résultat
    window_creation(m_window_T, "composition", m_width_T, m_height_T);

    shader_program(m_NoiseShaderProgram, TEMPO_PATH+"shaders/noise_shader.frag", TEMPO_PATH+"shaders/identity_shader.vert");
    display_quad(m_NoiseVAO);

    shader_program(m_ColorShaderProgram, TEMPO_PATH+"shaders/color_shader.frag", TEMPO_PATH+"shaders/identity_shader.vert");
    display_quad(m_ColorVAO);

    shader_program(m_CompositionShaderProgram, TEMPO_PATH+"shaders/composition_shader.frag", TEMPO_PATH+"shaders/identity_shader.vert");
    display_quad(m_CompositionVAO);


    // fbo noise 1
    if(!setFBO(m_fbo_N1, m_texture_N1, m_width_T, m_height_T)){
        return false;
    }
    drawNoise(m_fbo_N1, m_width_T, m_height_T, m_F1Princ, m_F1Spread, m_Or1Princ, m_Or1Spread, 1.);
    computeStatistiques(m_fbo_N1, m_width_T, m_height_T, m_mean_N1, m_var_N1);


    if(!setFBO(m_fbo_N1_ui, m_texture_N1_ui, m_width_N, m_height_N)){
        return false;
    }
    drawNoise(m_fbo_N1_ui, m_width_T, m_height_T, m_F1Princ, m_F1Spread, m_Or1Princ, m_Or1Spread, 1.);



    // fbo noise 2
    if(!setFBO(m_fbo_N2, m_texture_N2, m_width_T, m_height_T)){
        return false;
    }
    drawNoise(m_fbo_N2, m_width_T, m_height_T, m_F2Princ, m_F2Spread, m_Or2Princ, m_Or2Spread, 2.);
    computeStatistiques(m_fbo_N2, m_width_T, m_height_T, m_mean_N2, m_var_N2);

    if(!setFBO(m_fbo_N2_ui, m_texture_N2_ui, m_width_N, m_height_N)){
        return false;
    }
    drawNoise(m_fbo_N2_ui, m_width_T, m_height_T, m_F2Princ, m_F2Spread, m_Or2Princ, m_Or2Spread, 2.);



    // fbo color map
    if(!setFBO(m_fbo_H, m_texture_H, m_width_H, m_height_H)){
        return false;
    }

    // fbo composition
    if(!setFBO(m_fbo_T, m_texture_T, m_width_T, m_height_T)){
        return false;
    }
    computeProportions();

    initGui(m_window_T);



    return true;
}



void ccvt_application::onFinish() {
    terminateGui();

//    glDeleteVertexArrays(1, &m_PointsVAO);
//    glDeleteProgram(m_PointsShaderProgram);

    glDeleteVertexArrays(1, &m_ColorVAO);
    glDeleteProgram(m_ColorShaderProgram);
//
    glDeleteVertexArrays(1, &m_NoiseVAO);
    glDeleteProgram(m_NoiseShaderProgram);

    glDeleteVertexArrays(1, &m_CompositionVAO);
    glDeleteProgram(m_CompositionShaderProgram);

//    glfwDestroyWindow(m_window_H);
    glfwDestroyWindow(m_window_T);
    glfwTerminate();
}



bool ccvt_application::isRunning() {
    // close the window when "escape" is press
    bool is_escape_pressed = glfwGetKey(m_window_T, GLFW_KEY_ESCAPE) == GLFW_PRESS;
//                            or glfwGetKey(m_window_T, GLFW_KEY_ESCAPE) == GLFW_PRESS;

    if(is_escape_pressed){
//        glfwSetWindowShouldClose(m_window_H, true);
        glfwSetWindowShouldClose(m_window_T, true);
    }

    return !glfwWindowShouldClose(m_window_T);
}





void ccvt_application::onFrame() {

    ////////////////////////////////////////////////////////////////////
    glfwMakeContextCurrent(m_window_T);
    glViewport(0, 0, m_width_T, m_height_T); // (lower_left_x, lower_left_y, width, height)

    glfwPollEvents();



    // noise 1
    //////////////////////////////////////
    if(m_n1_changed){
        drawNoise(m_fbo_N1, m_width_T, m_height_T, m_F1Princ, m_F1Spread, m_Or1Princ, m_Or1Spread, 1.);
        drawNoise(m_fbo_N1_ui, m_width_T, m_height_T, m_F1Princ, m_F1Spread, m_Or1Princ, m_Or1Spread, 1.);
        m_n1_changed = false;
    }


    // noise 2
    //////////////////////////////////////
    if(m_n2_changed)
    {
        drawNoise(m_fbo_N2, m_width_T, m_height_T, m_F2Princ, m_F2Spread, m_Or2Princ, m_Or2Spread, 2.);
        drawNoise(m_fbo_N2_ui, m_width_T, m_height_T, m_F2Princ, m_F2Spread, m_Or2Princ, m_Or2Spread, 2.);
        m_n2_changed = false;
    }




    // color map
    //////////////////////////////////////
    drawCM();



    // composition
    //////////////////////////////////////
    drawComposition();

    updateGui();





    // swap the buffers
    glfwSwapBuffers(m_window_T);



}




////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ccvt_application::drawNoise(unsigned int &fbo, float width, float height, float freq_princ, float freq_spread, float omega_princ, float omega_spread, float seed) {
    glViewport(0, 0, width, height); // (lower_left_x, lower_left_y, width, height)
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    // clear color
    glClearColor(0.6f, 0.6f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // we're not using the stencil buffer now
    glEnable(GL_DEPTH_TEST);

    glUseProgram(m_NoiseShaderProgram);

    // unifoms
    glUniform2f(glGetUniformLocation(m_NoiseShaderProgram, "uRes"), width, height);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uFprinc"), freq_princ);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uFspread"), freq_spread);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uOprinc"), omega_princ);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uOspread"), omega_spread);
    glUniform1f(glGetUniformLocation(m_NoiseShaderProgram, "uSeed"), seed);

    glBindVertexArray(m_NoiseVAO);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}




void ccvt_application::drawCM() {
    glViewport(0, 0, m_width_H, m_height_H);
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_H);

    // clear color
    glClearColor(0.4f, 0.6f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // we're not using the stencil buffer now
    glEnable(GL_DEPTH_TEST);

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

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}





void ccvt_application::drawComposition() {
    glViewport(0, 0, m_width_T, m_height_T);
    glClearColor(0.9f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // we're not using the stencil buffer now
    glEnable(GL_DEPTH_TEST);

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
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ccvt_application::initGui(GLFWwindow* &window) {
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
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
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
}


void ccvt_application::cellPopup(int id)
{
    // identifiant de la cellule
    ImGui::Text(("cell "+std::to_string(id)).c_str());
    ImGui::SameLine(ImGui::GetContentRegionAvail().x-35);
    if (ImGui::Button("Close"))
    {
        ImGui::CloseCurrentPopup();
    }

    // position de la graine
    ImGui::Separator();
    ImGui::InputScalar("position x",   ImGuiDataType_Float,  &(m_points.at(id)._x),  NULL);
    ImGui::InputScalar("position y",   ImGuiDataType_Float,  &(m_points.at(id)._y),  NULL);

    // proportion
    ImVec4 current_color(m_cells.at(id)._r, m_cells.at(id)._g, m_cells.at(id)._b, 1.f);
    ImGui::PushStyleColor(ImGuiCol_FrameBg, current_color);
    ImGui::DragFloat("proportion", &m_cells.at(id)._cap, 0.005f,  0.01f, 1.f, "%f");
    if(ImGui::IsItemEdited()){normilizeCap();}
    ImGui::PopStyleColor(1);

    // couleur
    ImGui::Separator();
    static ImVec4 color;
    ImGui::ColorPicker4("ColorPicker", (float*)&color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel, &current_color.x);


    ImGui::GetStyle().ItemSpacing.x = 8.;
    if (ImGui::Button("Change color"))
    {
        m_cells.at(id)._r = color.x;
        m_cells.at(id)._g = color.y;
        m_cells.at(id)._b = color.z;
    }
    ImGui::SameLine();
    if (ImGui::Button("Duplicate"))
    {
        float posX = m_points.at(id)._x + 0.2*(float(std::rand()) / float(RAND_MAX)) - 0.1;
        float posY = m_points.at(id)._y + 0.2*(float(std::rand()) / float(RAND_MAX)) - 0.1;
        insertPoint(id, posX, posY,current_color.x, current_color.y, current_color.z);
    }
    ImGui::SameLine();
    if (ImGui::Button("Delete"))
    {
        deletePoint(id);
        ImGui::CloseCurrentPopup();
    }
    // ImGui::SameLine();
    // if (ImGui::Button("Close"))
    // {
    //     ImGui::CloseCurrentPopup();
    // }

    ImGui::EndPopup();
}



void ccvt_application::updateGui() {
    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    // ImGui::ShowDemoWindow(); // Show demo window! :)


    ImGui::Begin("Performances");
    ImGuiIO& io = ImGui::GetIO();
    ImGui::Text("Application average %.3f ms/frame", 1000.0f / io.Framerate);
    ImGui::Text("(%.3f FPS)", io.Framerate);
    ImGui::End();


    ImGui::Begin("Infos");
    if(ImGui::Button("Statistiques noises")){
        computeStatistiques(m_fbo_N1, m_width_T, m_height_T, m_mean_N1, m_var_N1);
        computeStatistiques(m_fbo_N2, m_width_T, m_height_T, m_mean_N2, m_var_N2);
    }
    ImGui::Text("Noise 1:");
    ImGui::BulletText(
            "Mean = %f\n",
            m_mean_N1
    );
    ImGui::BulletText(
            "Variance = %f\n",
            m_var_N1
    );

    ImGui::Text("Noise 2:");
    ImGui::BulletText(
            "Mean = %f\n",
            m_mean_N2
    );
    ImGui::BulletText(
            "Variance = %f\n",
            m_var_N2
    );

    ImGui::Separator();
    if(ImGui::Button("Proportion colors")){
        computeProportions();
    }
    // for(auto c=m_cells.begin(); c<m_cells.end(); c++){
    //     ImGui::ColorButton("color", ImVec4((*c)._r, (*c)._g, (*c)._b, 1.), ImGuiColorEditFlags_NoBorder);
    //     ImGui::SameLine();
    //     ImGui::Text(" %.2f%% (obj: %.2f%%) ", 100.*m_histo.at(c-m_cells.begin())._count/float(m_width_T*m_height_T), 100.*(*c)._cap);
    // }
    for(auto c=m_histo.begin(); c<m_histo.end(); c++){
        ImGui::ColorButton("color", ImVec4((*c)._r, (*c)._g, (*c)._b, 1.), ImGuiColorEditFlags_NoBorder);
        ImGui::SameLine();
        ImGui::Text(" %.2f%% (obj: %.2f%%) ", 100.*(*c)._count/float(m_width_T*m_height_T), 100.*(*c)._obj);
    }

    ImGui::Separator();
    static bool color_map = true;
    static bool noise_1 = true;
    static bool noise_2 = true;
    static bool composition = true;
    if(ImGui::Button("Save")){
        if(color_map){
            saveTexture(m_fbo_H, m_width_H, m_height_H, "color_map.png");
            saveData("color_map.txt");
        }
        if(noise_1){
            saveTexture(m_fbo_N1, m_width_T, m_height_T, "noise_1.png");
        }
        if(noise_2){
            saveTexture(m_fbo_N2, m_width_T, m_height_T, "noise_2.png");
        }
        if(composition){
            glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_T);
            drawComposition();
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            saveTexture(m_fbo_T, m_width_T, m_height_T, "composition.png");
        }
    }

    ImGui::Checkbox("color_map.ppm", &color_map);
    ImGui::Checkbox("noise_1.ppm", &noise_1);
    ImGui::Checkbox("noise_2.ppm", &noise_2);
    ImGui::Checkbox("composition.ppm", &composition);

    ImGui::Separator();
    ImGui::Text(m_infoBuffer.c_str());
    ImGui::End();






    ImGui::Begin("Color map");

    if(ImPlot::BeginPlot("map plot", ImVec2(m_width_H, m_height_H), ImPlotFlags_Equal | ImPlotFlags_NoLegend | ImPlotFlags_NoTitle))
    {
        ImPlot::SetupAxesLimits(0, 1, 0, 1, ImGuiCond_Always);
        ImPlot::SetupAxes(NULL,NULL,ImPlotAxisFlags_NoDecorations,ImPlotAxisFlags_NoDecorations);

        ImPlot::PlotImage("image", (void *)m_texture_H, ImPlotPoint(0, 0), ImPlotPoint(1, 1), ImVec2(0, 1), ImVec2(1, 0));


        // ajout de points
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(0) && ImGui::GetIO().KeyShift) { //  && ImGui::GetIO().KeyShift
            ImPlotPoint pt = ImPlot::GetPlotMousePos();
            addPoint(pt.x, pt.y, (float(std::rand()) / float(RAND_MAX)), (float(std::rand()) / float(RAND_MAX)), (float(std::rand()) / float(RAND_MAX)));
        }


        // déplacement
        for(auto p=m_points.begin(); p<m_points.end(); ++p){
            bool clicked;
            bool hovered;
            bool held;
            ImVec4 color = m_selected.active && m_selected.id==p-m_points.begin() ? ImVec4(1., 1., 1., 1.) : ImVec4(0., 0., 0., 1.);
            ImPlot::DragPoint(p-m_points.begin(), &(*p)._x, &(*p)._y, color, 6., ImPlotDragToolFlags_None, &clicked, &hovered, &held);
            if(hovered)
            {
                ImGui::BeginTooltip();
                ImGui::Text(("cell id : "+std::to_string(p-m_points.begin())).c_str());
                ImGui::Text(("position : "+std::to_string((*p)._x) +", " + std::to_string((*p)._y)).c_str());
                ImGui::EndTooltip();
            }
            if(clicked && ImGui::GetIO().KeyCtrl)
            {
                m_selected.id = p-m_points.begin();
                ImGui::OpenPopup("cell_param");
            }
            if(held)
            {
                m_selected.toogle = true;
            }
        }

        // graphe adjacence
        if(m_showAdjGraph){
            if(m_adjGraph_id.empty())
            {
                getCCVTadjGraph();
            }
            if(m_selected.toogle)
            {
                getCCVTadjGraph();
                m_selected.toogle = false;
            }
            for(int id=0; id<m_adjGraph_id.size()-1; id+=2)
            {
                int point_id = m_adjGraph_id.at(id);
                int point_id_n = m_adjGraph_id.at(id+1);
                std::vector<float> x{m_points.at(point_id)._x, m_points.at(point_id_n)._x};
                std::vector<float> y{m_points.at(point_id)._y, m_points.at(point_id_n)._y};
                ImPlot::SetNextLineStyle(ImVec4(1,0,0,1));
                ImPlot::PlotLine("adjacence", x.data(), y.data(), 2, 0, 0);
            }
        }



        // suppression et changement de couleur
        if (ImGui::BeginPopup("cell_param"))
        {
            int id = m_selected.id;
            cellPopup(id);
        }


        ImPlot::EndPlot();

    }






    // barre proportion des couleurs
    float window_x = ImGui::GetContentRegionAvail().x-m_CurrentPointsNb;
    ImGui::GetStyle().ItemSpacing.x = 1.;
    m_selected.active = false;
    for(auto c=m_cells.begin(); c<m_cells.end(); c++)
    {
        int id = c-m_cells.begin();

        if (c>m_cells.begin()){
            ImGui::SameLine();
        }

        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4((*c)._r, (*c)._g, (*c)._b, 1.f));
        if(ImGui::Button(("cell "+std::to_string(id)).c_str(), ImVec2(window_x*(*c)._cap, 60)))
        {
            m_selected.id = id;
            ImGui::OpenPopup("cell_param");

        }
        if(ImGui::IsItemHovered()){
            m_selected.id = id;
            m_selected.active = true;
        }
        ImGui::PopStyleColor(1);
    }




    // suppression et changement de couleur
    if (ImGui::BeginPopup("cell_param"))
    {
        int id = m_selected.id;
        cellPopup(id);
    }


    ImGui::Separator();
    ImGui::GetStyle().ItemSpacing.x = 8.;
    if(ImGui::Button("Same proportions")){
        equalizeCap();
    }
    ImGui::SameLine();
    if(ImGui::Button("Optimize")){
        optimizeCCVT();
    }
    ImGui::SameLine();
    ImGui::Checkbox("Show adjacence graphe", &m_showAdjGraph);

    ImGui::Separator();
    ImGui::InputScalar("size H",   ImGuiDataType_U16,  &m_size_density,  NULL);



    ImGui::Separator();
    ImGui::Text("Number of points :");
    ImGui::BulletText(
            "Max number = %i\n",
            m_MaxPointsNb
    );
    ImGui::BulletText(
            "Min number = %i\n",
            3
    );
    ImGui::BulletText(
            "Current number = %i\n",
            m_CurrentPointsNb
    );


    ImGui::End();










    ImGui::Begin("noise 1");

    ImGui::SliderFloat("F_0", &m_F1Princ, 10., 60.);
    if(ImGui::IsItemEdited()){m_n1_changed = true;}
    ImGui::SliderFloat("F spread", &m_F1Spread, 0., 20.);
    if(ImGui::IsItemEdited()){m_n1_changed = true;}
    ImGui::SliderAngle("O_0", &m_Or1Princ, 0., 180.);
    if(ImGui::IsItemEdited()){m_n1_changed = true;}
    ImGui::SliderAngle("O spread", &m_Or1Spread, 0., 90.);
    if(ImGui::IsItemEdited()){m_n1_changed = true;}


    // we access the ImGui window size
    const float window_width_N1 = m_width_N;//ImGui::GetContentRegionAvail().x;
    const float window_height_N1 = m_height_N;//ImGui::GetContentRegionAvail().y;

//    glViewport(0, 0, window_width_N1, window_height_N1);

    // we get the screen position of the window
    ImVec2 pos_N1 = ImGui::GetCursorScreenPos();

    // and here we can add our created texture as image to ImGui
    // unfortunately we need to use the cast to void* or I didn't find another way tbh
    ImGui::GetWindowDrawList()->AddImage(
            (void *)m_texture_N1_ui,
            ImVec2(pos_N1.x, pos_N1.y),
            ImVec2(pos_N1.x + window_width_N1, pos_N1.y + window_height_N1),
            ImVec2(0, 1),
            ImVec2(1, 0)
    );
    ImGui::End();






    ImGui::Begin("noise 2");

    ImGui::SliderFloat("F_0", &m_F2Princ, 10., 60.);
    if(ImGui::IsItemEdited()){m_n2_changed = true;}
    ImGui::SliderFloat("F spread", &m_F2Spread, 0., 20.);
    if(ImGui::IsItemEdited()){m_n2_changed = true;}
    ImGui::SliderAngle("O_0", &m_Or2Princ, 0., 180.);
    if(ImGui::IsItemEdited()){m_n2_changed = true;}
    ImGui::SliderAngle("O spread", &m_Or2Spread, 0., 90.);
    if(ImGui::IsItemEdited()){m_n2_changed = true;}


    // we access the ImGui window size
    const float window_width_N2 = m_width_N;//ImGui::GetContentRegionAvail().x;
    const float window_height_N2 = m_height_N;//ImGui::GetContentRegionAvail().y;

//    glViewport(0, 0, window_width_N2, window_height_N2);

    // we get the screen position of the window
    ImVec2 pos_N2 = ImGui::GetCursorScreenPos();

    // and here we can add our created texture as image to ImGui
    // unfortunately we need to use the cast to void* or I didn't find another way tbh
    ImGui::GetWindowDrawList()->AddImage(
            (void *)m_texture_N2_ui,
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
void ccvt_application::addPoint(float xPos, float yPos, float r, float g, float b) {
//    std::cout<<"new point : "<<xPos<<", "<<yPos<<std::endl;
    if(m_CurrentPointsNb < m_MaxPointsNb){
        float new_cap = 1.f/(m_cells.size()+1);

        m_points.emplace_back(xPos, yPos, 0.);
        m_cells.emplace_back(new_cap, r, g, b);
        m_CurrentPointsNb = m_points.size();

        normilizeCap();
        getCCVTadjGraph();
        // m_histo.resize(m_CurrentPointsNb, m_cells.end()->_cap);
        computeProportions();
    }

}

void ccvt_application::insertPoint(int vecPos, float xPos, float yPos, float r, float g, float b)
{
    if(m_CurrentPointsNb < m_MaxPointsNb)
    {
        float new_cap = m_cells.at(vecPos)._cap;// 1.f/(m_cells.size()+1);

        m_points.insert(m_points.begin()+vecPos, point_info{xPos, yPos, 0.}); //emplace_back(xPos, yPos, 0.);
        m_cells.insert(m_cells.begin()+vecPos, cell_info{new_cap, r, g, b}); //emplace_back(new_cap, r, g, b);
        m_CurrentPointsNb = m_points.size();


        normilizeCap();
        getCCVTadjGraph();
        // m_histo.resize(m_CurrentPointsNb, m_cells.end()->_cap);
        computeProportions();
    }
}


void ccvt_application::deletePoint(int id) {
    if(m_CurrentPointsNb>3){
        m_points.erase(m_points.begin()+id);
        m_cells.erase(m_cells.begin()+id);
        m_histo.erase(m_histo.begin()+id);
        m_CurrentPointsNb = m_points.size();

        normilizeCap();
        getCCVTadjGraph();
        computeProportions();
    }
}

auto sum_caps = [](float a, cell_info b){return a + b._cap;};
void ccvt_application::normilizeCap() {
    float sum = std::accumulate(m_cells.begin(), m_cells.end(), 0.f, sum_caps);
    for(auto &ci:m_cells){
        ci._cap = ci._cap/sum;
    }
}

void ccvt_application::equalizeCap() {
    for(auto &ci:m_cells){
        ci._cap = 1./m_CurrentPointsNb;
    }
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


void ccvt_application::getCCVTadjGraph()
{
    updateCCVT();
    m_ccvt.get_adjacence_graph(m_adjGraph_id);
}





void ccvt_application::optimizeCCVT(){
    normilizeCap();
    updateCCVT();
    std::cout<<"optimizing CCVT..."<<std::endl;
    FT stepX = 0.05; // pour les positions
    FT stepW = 0.1; // pour les poids
    FT epsilon = 1.;
    unsigned max_newton_iters = 500;
    unsigned max_iters = 50;

    // int iter_opt = m_ccvt.optimize_H(stepW, stepX, max_newton_iters, epsilon, max_iters);
    Timer::start_timer(m_timer, COLOR_BLUE, "optimize proportions");
    unsigned int iter_opt = m_ccvt.optimize_all(stepW, stepX, max_newton_iters, epsilon, max_iters, std::cout);
    double duration = Timer::stop_timer(m_timer, COLOR_BLUE);

    std::cout<<iter_opt<<" itérations (max : "<<max_newton_iters*max_iters*2.<<")"<<std::endl;
    //m_infoBuffer = std::to_string(iter_opt) + " iterations in "+std::to_string(duration) + "s (max : " + std::to_string(max_newton_iters) + ")\n";
    std::sprintf(m_infoBuffer.data(), "%i iteration done in %.3fs", iter_opt, duration);

    getCCVTcells();
    getCCVTadjGraph();
    computeProportions();
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

    computeStatistiques(m_fbo_N1, m_width_T, m_height_T, m_mean_N1, m_var_N1);
    computeStatistiques(m_fbo_N2, m_width_T, m_height_T, m_mean_N2, m_var_N2);

    m_ccvt.set_domain(m_mean_N1, m_mean_N2, m_var_N1, m_var_N2, m_size_density, m_size_density, m_max_val);
    m_ccvt.set_custom_proportions(caps);
    m_ccvt.set_sites(points, weights);
    m_ccvt.set_colors(R, G, B);

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ccvt_application::saveTexture(unsigned int fbo_id, int width, int height, const std::string &filename) {
    std::cout<<"saving "<<filename<<"..."<<std::endl;

    glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);

    glViewport(0, 0, width, height);

    std::vector<float> pixel(3*width*height);
    glReadPixels(0, 0, width, height, GL_RGB, GL_FLOAT, pixel.data());


    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    std::string fn = TEMPO_PATH+"results/CCVT_CGAL/application/"+filename; 

	stbi_flip_vertically_on_write(0);
	stbi_write_png(fn.c_str(), width, height, 3, pixel.data(), width * 3);

}


void ccvt_application::saveData(const std::string &filename)
{
    std::string fn = TEMPO_PATH+"results/CCVT_CGAL/application/"+filename;  
    std::ofstream output(fn.c_str());
    for(auto p=m_points.begin(); p<m_points.end(); p++)
    {
        output<<"position: ("<<(*p)._x<<", "<<(*p)._y<<"), ";
        output<<"weigth: "<<(*p)._w<<", ";
        output<<"colour: ("<<m_cells.at(p-m_points.begin())._r<<", "<<m_cells.at(p-m_points.begin())._g<<", "<<m_cells.at(p-m_points.begin())._b<<").";
        output<<std::endl;
    }
    output.close();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ccvt_application::computeStatistiques(unsigned int fbo_id, int width, int height, float &mean, float &var) {
    glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);

    glViewport(0, 0, width, height);

    std::vector<float> pixel(width*height);
    glReadPixels(0, 0, width, height, GL_RED, GL_FLOAT, pixel.data());

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    float moy = computeMean(pixel);
    float moyS = computeSquareMean(pixel);

    mean = moy;
    var = moyS - moy*moy;
}


float ccvt_application::computeMean(const std::vector<float>& data) {
    float value = 0;
    for(auto d:data){
        value += d;
    }
    return value/data.size();
}


float ccvt_application::computeSquareMean(const std::vector<float>& data) {
    float value = 0;
    for(auto d:data){
        value += d*d;
    }
    return value/data.size();
}

void ccvt_application::computeProportions() {
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_T);
    drawComposition();

    std::vector<float> pixel(3*m_width_T*m_height_T);
    glReadPixels(0, 0, m_width_T, m_height_T, GL_RGB, GL_FLOAT, pixel.data());

    glBindFramebuffer(GL_FRAMEBUFFER, 0);


    std::vector<histo_info> histo;//(m_CurrentPointsNb);
    for(int p=0; p<m_width_T*m_height_T*3; p+=3)
    {
        bool found = false;
        for(auto c = histo.begin(); c<histo.end(); ++c)
        {
            found = abs(pixel[p]-(*c)._r)<0.01 && abs(pixel[p+1]-(*c)._g)<0.01 && abs(pixel[p+2]-(*c)._b)<0.01;
            if(found)
            {
                histo.at(c-histo.begin()).incr();
                break;
            }
        }
        if(!found)
        {
            histo.emplace_back(pixel[p], pixel[p+1], pixel[p+2]);
        }


        // for(auto c=m_cells.begin(); c<m_cells.end(); c++){
        //
        //     bool found = abs(pixel[p]-(*c)._r)<0.01 and abs(pixel[p+1]-(*c)._g)<0.01 and abs(pixel[p+2]-(*c)._b)<0.01;
        //     if(found){
        //         histo.at(c-m_cells.begin()) += 1;
        //     }
        // }
    }

    for(auto c=m_cells.begin(); c<m_cells.end(); ++c)
    {
        for(auto h = histo.begin(); h<histo.end(); ++h)
        {
            bool found = abs((*h)._r-(*c)._r)<0.01 && abs((*h)._g-(*c)._g)<0.01 && abs((*h)._b-(*c)._b)<0.01;
            if(found)
            {
                (*h)._obj+=(*c)._cap;
            }
        }
    }


    m_histo=histo;
}