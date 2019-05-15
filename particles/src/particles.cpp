#include "particles.h"

#include <math.h>
#include "gl/shaders/ShaderAttribLocations.h"

Particles::Particles() : m_numRows(100), m_numCols(m_numRows)
{
}


/**
 * Returns a pseudo-random value between -1.0 and 1.0 for the given row and column.
 */
float Particles::randValue(int row, int col) {
    return -1.0 + 2.0 * glm::fract(sin(row * 127.1f + col * 311.7f) * 43758.5453123f);
}


/**
 * Returns the object-space position for the Particles vertex at the given row and column.
 */
glm::vec3 Particles::getPosition(int row, int col) {
    glm::vec3 position;

    position.x = 10 * row/m_numRows - 5;
    position.y = getHeight(row, col, 20);
    position.z = 10 * col/m_numCols - 5;

    position.y += 0.4*getHeight(row, col, 10);
    position.y += 0.3*getHeight(row, col, 4);

    return position;
}

float Particles::getHeight(int row, int col, int frequency) {
    float A = randValue(row / frequency, col / frequency);
    float B = randValue(row / frequency, col / frequency + 1);
    float C = randValue(row / frequency + 1, col / frequency);
    float D = randValue(row / frequency + 1, col / frequency + 1);

    float upDist = glm::fract(row / static_cast<float>(frequency));
    float leftDist = glm::fract(col / static_cast<float>(frequency));
    upDist = 3 * std::pow(upDist, 2) - 2 * std::pow(upDist, 3);
    leftDist = 3 * std::pow(leftDist, 2) - 2 * std::pow(leftDist, 3);

    float interpHeight = (1.0 - upDist) * ((1.0 - leftDist) * A + leftDist * B) + upDist * ((1.0 - leftDist) * C + leftDist * D);
    return interpHeight;
}


/**
 * Returns the normal vector for the Particles vertex at the given row and column.
 */
glm::vec3 Particles::getNormal(int row, int col) {
    // TODO: Compute the normal at the given row and column using the positions of the
    //       neighboring vertices.
    glm::vec3 p = glm::vec3(getPosition(row, col));
    glm::vec3 n0 = glm::vec3(getPosition(row, col + 1) - p);
    glm::vec3 n1 = glm::vec3(getPosition(row - 1, col + 1) - p);
    glm::vec3 n2 = glm::vec3(getPosition(row - 1, col) - p);
    glm::vec3 n3 = glm::vec3(getPosition(row - 1, col - 1) - p);
    glm::vec3 n4 = glm::vec3(getPosition(row, col - 1) - p);
    glm::vec3 n5 = glm::vec3(getPosition(row + 1, col - 1) - p);
    glm::vec3 n6 = glm::vec3(getPosition(row + 1, col) - p);
    glm::vec3 n7 = glm::vec3(getPosition(row + 1, col + 1) - p);

    glm::vec3 norm1 = glm::normalize(glm::cross(n1, n0));
    glm::vec3 norm2 = glm::normalize(glm::cross(n2, n1));
    glm::vec3 norm3 = glm::normalize(glm::cross(n3, n2));
    glm::vec3 norm4 = glm::normalize(glm::cross(n4, n3));
    glm::vec3 norm5 = glm::normalize(glm::cross(n5, n4));
    glm::vec3 norm6 = glm::normalize(glm::cross(n6, n5));
    glm::vec3 norm7 = glm::normalize(glm::cross(n7, n6));
    glm::vec3 norm8 = glm::normalize(glm::cross(n0, n7));

    glm::vec3 average = (norm1 + norm2 + norm3 + norm4 + norm5 + norm6 + norm7 + norm8) / 8.f;

    return average;
}


/**
 * Initializes the Particles by storing positions and normals in a vertex buffer.
 */
void Particles::init() {
    // TODO: Change from GL_LINE to GL_FILL in order to render full triangles instead of wireframe.
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


    // Initializes a grid of vertices using triangle strips.
    int numVertices = (m_numRows - 1) * (2 * m_numCols + 2);
    std::vector<glm::vec3> data(2 * numVertices);
    int index = 0;
    for (int row = 0; row < m_numRows - 1; row++) {
        for (int col = m_numCols - 1; col >= 0; col--) {
            data[index++] = getPosition(row, col);
            data[index++] = getNormal  (row, col);
            data[index++] = getPosition(row + 1, col);
            data[index++] = getNormal  (row + 1, col);
        }
        data[index++] = getPosition(row + 1, 0);
        data[index++] = getNormal  (row + 1, 0);
        data[index++] = getPosition(row + 1, m_numCols - 1);
        data[index++] = getNormal  (row + 1, m_numCols - 1);
    }

    // Initialize OpenGLShape.
    m_shape = std::make_unique<OpenGLShape>();
    m_shape->setVertexData(&data[0][0], data.size() * 3, VBO::GEOMETRY_LAYOUT::LAYOUT_TRIANGLE_STRIP, numVertices);
    m_shape->setAttribute(ShaderAttrib::POSITION, 3, 0, VBOAttribMarker::DATA_TYPE::FLOAT, false);
    m_shape->setAttribute(ShaderAttrib::NORMAL, 3, sizeof(glm::vec3), VBOAttribMarker::DATA_TYPE::FLOAT, false);
    m_shape->buildVAO();
}


/**
 * Draws the Particles.
 */
void Particles::draw()
{
    m_shape->draw();
}
