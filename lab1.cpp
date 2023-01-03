#include <random>

#include <spdlog/spdlog.h>

#include "Labs/1-Drawing2D/tasks.h"

#include<iostream>
using namespace std;

using VCX::Labs::Common::ImageRGB;

namespace VCX::Labs::Drawing2D {
    /******************* 1.Image Dithering *****************/
    void DitheringThreshold(
        ImageRGB &       output,
        ImageRGB const & input) {
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input[{ x, y }];
                output.SetAt({ x, y }, {
                                           color.r > 0.5 ? 1 : 0,
                                           color.g > 0.5 ? 1 : 0,
                                           color.b > 0.5 ? 1 : 0,
                                       });
            }
    }

    void DitheringRandomUniform(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                float noise = rand() / (float)RAND_MAX - 0.5;
                glm::vec3 color = input[{ x, y }];
                output.SetAt({ x, y }, {
                                           (color.r + noise > 0.5) ? 1 : 0,
                                           (color.g + noise > 0.5) ? 1 : 0,
                                           (color.b + noise > 0.5) ? 1 : 0,
                                       });
            }
    }

    void DitheringRandomBlueNoise(
        ImageRGB &       output,
        ImageRGB const & input,
        ImageRGB const & noise) {
        // your code here:
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                float r_offset = rand() / (float)RAND_MAX - 0.5;
                glm::vec3 color = input[{ x, y }];
                glm::vec3 color_noise = noise[{ x, y }];
                output.SetAt({ x, y }, {
                                           color.r + color_noise.r - 0.5,
                                           color.g + color_noise.g - 0.5,
                                           color.b + color_noise.b - 0.5,
                                       });
            }
    }

    void DitheringOrdered(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        float matrix[3][3] = {{0.6, 0.8, 0.4}, {0.1, 0, 0.3}, {0.5, 0.2, 0.7}};
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input[{ x, y }];
                for(std::size_t dx = 0; dx < 3; ++dx) 
                    for(std::size_t dy = 0; dy < 3; ++dy) {
                        output.SetAt({ 3 * x + dx, 3 * y + dy }, {
                                           color.r > matrix[dx][dy] ? 1: 0,
                                           color.g > matrix[dx][dy] ? 1: 0,
                                           color.b > matrix[dx][dy] ? 1: 0,
                                       });   
                    }
            }
    }

    void DitheringErrorDiffuse(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        ImageRGB tmp = Common::CreatePureImageRGB(input.GetSizeX(), input.GetSizeY(), {0,0,0});
        int dx[4] = {1, -1, 0, 1};
        int dy[4] = {0, 1, 1, 1};
        float ratio[4] = {7.0/16.0, 3.0/16.0, 5.0/16, 1.0/16.0};
        for (std::size_t y = 0; y < input.GetSizeY(); ++y)
            for (std::size_t x = 0; x < input.GetSizeX(); ++x){
                glm::vec3 color = input[{x, y}];
                float color_r = color.r + tmp[{x, y}].r, color_g = color.g + tmp[{x, y}].g, color_b = color.b + tmp[{x, y}].b;
                output.SetAt({ x, y }, {
                    color_r > 0.5 ? 1 : 0,
                    color_g > 0.5 ? 1 : 0,
                    color_b > 0.5 ? 1 : 0,
                });
                for(int i = 0; i < 4; i++) {
                    int xx = x + dx[i];
                    int yy = y + dy[i];
                    if(xx >= 0 && yy >= 0 && xx < input.GetSizeX() && yy < input.GetSizeY()) {
                        glm::vec3 tmp_color = tmp[{ (std::size_t)xx, (std::size_t)yy }];
                        tmp.SetAt({ (std::size_t)xx, (std::size_t)yy }, {
                            tmp_color.r + (color_r - (color_r > 0.5 ? 1 : 0)) * ratio[i],
                            tmp_color.g + (color_g - (color_g > 0.5 ? 1 : 0)) * ratio[i],
                            tmp_color.b + (color_b - (color_b > 0.5 ? 1 : 0)) * ratio[i],
                        });
                    }
                }
            }
    }

    /******************* 2.Image Filtering *****************/
    void Blur(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        int dx[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
        int dy[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                float color_r = 0.0, color_g = 0.0, color_b = 0.0;
                int cnt = 0;
                for(int i = 0; i < 9; i++) {
                    int xx = x + dx[i], yy = y + dy[i];
                    if(xx >= 0 && yy >= 0 && xx < input.GetSizeX() && yy < input.GetSizeY()) {
                        cnt++;
                        glm::vec3 tmp_color = input[{(std::size_t)xx, (std::size_t)yy}];
                        color_r += tmp_color.r;
                        color_g += tmp_color.g;
                        color_b += tmp_color.b;
                    }   
                }
                output.SetAt({x, y}, {
                    color_r / cnt,
                    color_g / cnt,
                    color_b / cnt,
                });
            }
    }

    void Edge(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        int dx[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
        int dy[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
        int sobel_x[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
        int sobel_y[9] = {1, 0, -1, 2, 0, -2, 1, 0, -1};
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                float color[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
                for(int i = 0; i < 9; i++) {
                    int xx = x + dx[i], yy = y + dy[i];
                    if(xx >= 0 && yy >= 0 && xx < input.GetSizeX() && yy < input.GetSizeY()) {
                        color[0][0] += input[{(std::size_t)xx, (std::size_t)yy}].r * sobel_y[i];
                        color[0][1] += input[{(std::size_t)xx, (std::size_t)yy}].g * sobel_y[i];
                        color[0][2] += input[{(std::size_t)xx, (std::size_t)yy}].b * sobel_y[i];
                        color[1][0] += input[{(std::size_t)xx, (std::size_t)yy}].r * sobel_x[i];
                        color[1][1] += input[{(std::size_t)xx, (std::size_t)yy}].g * sobel_x[i];
                        color[1][2] += input[{(std::size_t)xx, (std::size_t)yy}].b * sobel_x[i];
                    }   
                }
                output.SetAt({x, y}, {
                    sqrt(color[0][0] * color[0][0] + color[1][0] * color[1][0]) / 2,
                    sqrt(color[0][1] * color[0][1] + color[1][1] * color[1][1]) / 2,
                    sqrt(color[0][2] * color[0][2] + color[1][2] * color[1][2]) / 2,
                });

            }
    }

    /******************* 3. Image Inpainting *****************/
    void Inpainting(
        ImageRGB &         output,
        ImageRGB const &   inputBack,
        ImageRGB const &   inputFront,
        const glm::ivec2 & offset) {
        output             = inputBack;
        size_t      width  = inputFront.GetSizeX();
        size_t      height = inputFront.GetSizeY();
        glm::vec3 * g      = new glm::vec3[width * height];
        memset(g, 0, sizeof(glm::vec3) * width * height);
        // set boundary condition
        for (std::size_t y = 0; y < height; ++y) {
            // set boundary for (0, y), your code: g[y * width] = ?
            g[y * width] = inputBack[{(std::size_t)offset.x, y + offset.y}] - inputFront[{0, y}];
            // set boundary for (width - 1, y), your code: g[y * width + width - 1] = ?
            g[y * width + width - 1] = inputBack[{width - 1 + offset.x, y + offset.y}] - inputFront[{width - 1, y}];

        }
        for (std::size_t x = 0; x < width; ++x) {
            // set boundary for (x, 0), your code: g[x] = ?
            g[x] = inputBack[{x + offset.x, (std::size_t)offset.y}] - inputFront[{x, 0}];
            // set boundary for (x, height - 1), your code: g[(height - 1) * width + x] = ?
            g[(height - 1) * width + x] = inputBack[{x + offset.x, height - 1 + offset.y}] - inputFront[{x, height - 1}];
        }

        // Jacobi iteration, solve Ag = b
        for (int iter = 0; iter < 8000; ++iter) {
            for (std::size_t y = 1; y < height - 1; ++y)
                for (std::size_t x = 1; x < width - 1; ++x) {
                    g[y * width + x] = (g[(y - 1) * width + x] + g[(y + 1) * width + x] + g[y * width + x - 1] + g[y * width + x + 1]);
                    g[y * width + x] = g[y * width + x] * glm::vec3(0.25);
                }
        }

        for (std::size_t y = 0; y < inputFront.GetSizeY(); ++y)
            for (std::size_t x = 0; x < inputFront.GetSizeX(); ++x) {
                glm::vec3 color = g[y * width + x] + inputFront.GetAt({ x, y });
                output.SetAt({ x + offset.x, y + offset.y }, color);
            }
        delete[] g;
    }

    /******************* 4. Line Drawing *****************/
    void DrawLine(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1) {
        // your code here:
        int x0 = p0.x < p1.x? p0.x: p1.x;
        int y0 = p0.y < p1.y? p0.y: p1.y;
        int x1 = p0.x < p1.x? p1.x: p0.x;
        int y1 = p0.y < p1.y? p1.y: p0.y;
        int ddx = 2*(x1-x0), ddy = 2*(y1-y0);
        //cout << ddx << ";" << ddy << endl;
        if (ddx <= 1){
            if (y0 > y1){
                int tmp = y0; y0 = y1; y1 = tmp;
            }
            for (int y = y0; y <= y1; ++y)
                canvas.SetAt({(std::size_t)x0, (std::size_t)y}, color);
            return ;
        }
        float m = abs((float)ddy / ddx);
        if (ddy >= 0 && m < 1){
            int F = ddy - ddx / 2;
            for (int x = x0, y = y0; x <= x1; ++x){
                F += ddy;
                if (F >= 0){
                    F -= ddx;
                    y++;
                }
                canvas.SetAt({(std::size_t)x, (std::size_t)y}, color);
            }
        }
        else if (ddy >= 0 && m >= 1){
            int F = ddy / 2 - ddx;
            for (int x = x0, y = y0; x<= x1; ++y){
                F -= ddx;
                if (F < 0){
                    F += ddy;
                    x++;
                }
                canvas.SetAt({(std::size_t)x, (std::size_t)y}, color);
            }
        }
        else if (ddy < 0 && m < 1){
            int F = ddy + ddx / 2;
            for (int x = x0, y = y0; x <= x1; ++x){
                F += ddy;
                if (F < 0){
                    F += ddx;
                    y--;
                }
                canvas.SetAt({(std::size_t)x, (std::size_t)y}, color);
            }
        }
        else{
            int F = ddy / 2 + ddx;
            for (int x = x0, y = y0; x <= x1; --y){
                F += ddx;
                if (F >= 0){
                    F += ddy;
                    x++;
                }
                canvas.SetAt({(std::size_t)x, (std::size_t)y}, color);
            }
        }
    }

    /******************* 5. Triangle Drawing *****************/
    void swap(int & x, int & y) {
        int t = x;
        x = y;
        y = t;
    }
    void DrawTriangleFilled(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1,
        glm::ivec2 const p2) {
        // your code here:
        int x[3], y[3];
        x[0] = p0.x; y[0] = p0.y;
        x[1] = p1.x; y[1] = p1.y;
        x[2] = p2.x; y[2] = p2.y;
        // 按照y从小到大排列
        if(y[0] > y[1]) {
            swap(x[0], x[1]);
            swap(y[0], y[1]);
        }
        if(y[0] > y[2]) {
            swap(x[0], x[2]);
            swap(y[0], y[2]);
        }
        if(y[1] > y[2]) {
            swap(x[1], x[2]);
            swap(y[1], y[2]);
        }
        if(!(y[0] <= y[1] && y[1] <= y[2])) return;
        float tmp = (float)(y[1] - y[0]) / (y[2] - y[0]) * (x[2] - x[0]) + x[0];
        float dx1 = (float)(x[1]-x[0]) / (y[1]-y[0]);
        float dx2 = (float)(tmp-x[0]) / (y[1]-y[0]);
        float x1 = x[0], x2 = x[0];
        for(std::size_t yy = y[0]; yy <= y[1]; yy++) {
            for(std::size_t xx = (std::size_t)x1; xx <= (std::size_t)x2; xx++) {
                canvas.SetAt({xx, yy}, color);
            }
            x1 += dx1;
            x2 += dx2;
        }
        dx1 = (float)(x[1]-x[2]) / (y[1]-y[2]);
        dx2 = (float)(tmp-x[2]) / (y[1]-y[2]);
        x1 = x[2], x2 = x[2];
        for(std::size_t yy = y[2]; yy >= y[1]; yy--) {
            for(std::size_t xx = (std::size_t)x1; xx <= (std::size_t)x2; xx++) {
                canvas.SetAt({xx, yy}, color);
            }
            x1 -= dx1;
            x2 -= dx2;
        }
    }

    /******************* 6. Image Supersampling *****************/
    void Supersample(
        ImageRGB &       output,
        ImageRGB const & input,
        int              rate) {
        // your code here:
        std::size_t length = output.GetSizeX() * rate, height = output.GetSizeY() * rate;
        std::size_t length_org = input.GetSizeX(), height_org = input.GetSizeY();
        ImageRGB temp = Common::CreatePureImageRGB(length, height, { 0, 0, 0 });
        for(std::size_t x = 0; x < length; x++) 
            for(std::size_t y = 0; y < height; y++) {
                temp.SetAt({x, y}, input[{(std::size_t)((float)x / length * length_org), (std::size_t)((float)y / height * height_org)}]);
            }
        for(std::size_t x = 0; x < output.GetSizeX(); x++) 
            for(std::size_t y = 0; y < output.GetSizeY(); y++) {
                float color_r = 0.0, color_g = 0.0, color_b = 0.0;
                for(std::size_t xx = 0; xx < rate; xx++)
                    for(std::size_t yy = 0; yy < rate; yy++) {
                        color_r += temp[{x * rate + xx, y * rate + yy}].r;
                        color_g += temp[{x * rate + xx, y * rate + yy}].g;
                        color_b += temp[{x * rate + xx, y * rate + yy}].b;
                    }
                color_r /= (rate * rate);
                color_g /= (rate * rate);
                color_b /= (rate * rate);
                output.SetAt({x, y}, {color_r, color_g, color_b});
            }
    }

    /******************* 7. Bezier Curve *****************/
    glm::vec2 CalculateBezierPoint(
        std::span<glm::vec2> points,
        float const          t) {
        // your code here:
       int n = points.size();

        glm::vec2 * temp = new glm::vec2[n];
        for(int i = 0; i < n; i++) {
            temp[i].x = points[i].x;
            temp[i].y = points[i].y;
        }
        for(int k = n - 1; k > 0; k--)
            for(int i = 0; i < k; i++) {
                temp[i].x = (1 - t) * temp[i].x + t * temp[i + 1].x;
                temp[i].y = (1 - t) * temp[i].y + t * temp[i + 1].y;
            }
        return glm::vec2 {temp[0].x, temp[0].y};
    }
} // namespace VCX::Labs::Drawing2D