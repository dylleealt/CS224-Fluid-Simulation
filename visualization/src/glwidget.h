#ifndef GLWIDGET_H
#define GLWIDGET_H
#include "GL/glew.h"
#include <QGLWidget>
#include <QTimer>

#include "glm/glm.hpp"            // glm::vec*, mat*, and basic glm functions
#include "glm/gtx/transform.hpp"  // glm::translate, scale, rotate
#include "glm/gtc/type_ptr.hpp"   // glm::value_ptr

#include <memory>  // std::unique_ptr

#include "gl/datatype/FBO.h"

class OpenGLShape;

using namespace CS123::GL;

class GLWidget : public QGLWidget {
    Q_OBJECT

public:
    GLWidget(QGLFormat format, QWidget *parent = 0);
    ~GLWidget();

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);
    void mousePressEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *e);

private:
    void drawBlur();
    void drawParticles();
    void drawFluidParticles();
    void setParticleViewport();

    void rebuildMatrices();

    int m_width;
    int m_height;

    GLuint m_phongProgram;
    GLuint m_textureProgram;
    GLuint m_horizontalBlurProgram;
    GLuint m_verticalBlurProgram;
    GLuint m_particleUpdateProgram;
    GLuint m_particleDrawProgram;
    GLuint m_fluidsUpdateProgram;
    GLuint m_fluidsDrawProgram;

    std::unique_ptr<OpenGLShape> m_quad;
    std::unique_ptr<OpenGLShape> m_sphere;

    std::unique_ptr<FBO> m_blurFBO1;
    std::unique_ptr<FBO> m_blurFBO2;

    GLuint m_particlesVAO;
    std::shared_ptr<FBO> m_particlesFBO1;
    std::shared_ptr<FBO> m_particlesFBO2;
    int m_numParticles;

    GLuint m_fluidsVAO;
    std::shared_ptr<FBO> m_fluidsFBO1;
    std::shared_ptr<FBO> m_fluidsFBO2;
    int m_sizeX;
    int m_sizeY;
    int m_sizeZ;
    int m_gridSize;
    bool m_firstPass;
    bool m_evenPass;

    glm::mat4 m_view, m_projection;

    /** For mouse interaction. */
    float m_angleX, m_angleY, m_zoom;
    QPoint m_prevMousePos;
};

#endif // GLWIDGET_H
