#include "glwidget.h"

#include "cs123_lib/resourceloader.h"
#include "cs123_lib/errorchecker.h"
#include <QMouseEvent>
#include <QWheelEvent>
#include <iostream>

#define PI 3.14159265f

GLWidget::GLWidget(QGLFormat format, QWidget *parent)
    : QGLWidget(format, parent), m_angleX(-0.9f), m_angleY(-2.4f), m_zoom(10.f)
{}

GLWidget::~GLWidget()
{}

void GLWidget::initializeGL() {
    ResourceLoader::initializeGlew();
    resizeGL(width(), height());

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    // Set the color to set the screen when the color buffer is cleared.
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    m_program = ResourceLoader::createShaderProgram(":/shaders/shader.vert", ":/shaders/shader.frag");
    m_particles.init();

    rebuildMatrices();
    paint();
}

void GLWidget::paint() {
//    // Bind shader program.
//    glUseProgram(m_program);

//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    // Set uniforms.
//    glUniformMatrix4fv(glGetUniformLocation(m_program, "model"), 1, GL_FALSE, glm::value_ptr(m_model));
//    glUniformMatrix4fv(glGetUniformLocation(m_program, "view"), 1, GL_FALSE, glm::value_ptr(m_view));
//    glUniformMatrix4fv(glGetUniformLocation(m_program, "projection"), 1, GL_FALSE, glm::value_ptr(m_projection));

//    // m_particles.draw();
//    this->grab().save("img.png");
//    std::cout<<"here"<<std::endl;
//    int frameNum = m_particles.getFrameNumber();
//    QPixmap pixmap(this->size());
//    this->render(&pixmap);
//    pixmap.save(QString("img").append(QString::number(frameNum)).append(QString(".png")));

//    // Unbind shader program.
//    glUseProgram(0);
}

void GLWidget::paintGL() {
    // Bind shader program.
    glUseProgram(m_program);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Set uniforms.
    glUniformMatrix4fv(glGetUniformLocation(m_program, "model"), 1, GL_FALSE, glm::value_ptr(m_model));
    glUniformMatrix4fv(glGetUniformLocation(m_program, "view"), 1, GL_FALSE, glm::value_ptr(m_view));
    glUniformMatrix4fv(glGetUniformLocation(m_program, "projection"), 1, GL_FALSE, glm::value_ptr(m_projection));

    m_particles.draw();
    int frameNumber = m_particles.getFrameNumber() + 1;
    glFlush();
    if (frameNumber < 10){
        grabFrameBuffer().save(QString("img/img000").append(QString::number(frameNumber)).append(QString(".png")));
    } else if (frameNumber < 100) {
        grabFrameBuffer().save(QString("img/img00").append(QString::number(frameNumber)).append(QString(".png")));
    } else if (frameNumber < 1000) {
        grabFrameBuffer().save(QString("img/img0").append(QString::number(frameNumber)).append(QString(".png")));
    }

    // Unbind shader program.
    glUseProgram(0);
}

void GLWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
    rebuildMatrices();
}

void GLWidget::mousePressEvent(QMouseEvent *event) {
    m_prevMousePos = event->pos();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
    m_angleX += 10 * (event->x() - m_prevMousePos.x()) / (float) width();
    m_angleY += 10 * (event->y() - m_prevMousePos.y()) / (float) height();
    m_prevMousePos = event->pos();
    rebuildMatrices();
}

void GLWidget::wheelEvent(QWheelEvent *event) {
    m_zoom -= event->delta() / 100.f;
    rebuildMatrices();
}

void GLWidget::rebuildMatrices() {
    m_model = glm::mat4(1.f);
    m_view = glm::translate(glm::vec3(0, 0, -3.f * m_zoom)) *
             glm::rotate(m_angleY, glm::vec3(1,0,0)) *
             glm::rotate(m_angleX, glm::vec3(0,1,0));
    m_projection = glm::perspective(0.8f, (float)width()/height(), 0.1f, 100.f);
    update();
}
