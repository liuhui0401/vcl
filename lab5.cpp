#include "Labs/5-Visualization/tasks.h"

#include <numbers>
#include <iostream>
#include <random>
#include <stdlib.h>
#include <algorithm>

using VCX::Labs::Common::ImageRGB;
namespace VCX::Labs::Visualization {

    struct CoordinateStates {
        // your code here
    };

    bool cmp(const Car& a, const Car& b) {
        return (a.mileage > b.mileage);
    }

    bool PaintParallelCoordinates(Common::ImageRGB & input, InteractProxy const & proxy, std::vector<Car> const & data, bool force) {
        // your code here
        // for example: 
        //   static CoordinateStates states(data);
        //   SetBackGround(input, glm::vec4(1));
        //   ...
        SetBackGround(input, glm::vec4(1));
       
        //cylinders
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(137.0 / 1920.0, 97.0 / 1920.0), 40.0 / 1920.0, "cylinders")) {
            std::cout << "wrong print" << '\n';
        }
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(137.0 / 1920.0, 155.0 / 1920.0), 40.0 / 1920.0, "9")) {
            std::cout << "wrong print" << '\n';
        }
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(118.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        DrawLine(input, glm::vec4(0.4, 0.4, 0.4, 0.7), glm::vec2(137.0 / 1920.0, 193.0 / 1920.0), glm::vec2(137.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(138.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(137.0 / 1920.0, 1764.0 / 1920.0), 40.0 / 1920.0, "2")) {
            std::cout << "wrong print" << '\n';
        }

        //displacement
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(411.0 / 1920.0, 97.0 / 1920.0), 40.0 / 1920.0, "displacement")) {
            std::cout << "wrong print" << '\n';
        }
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(411.0 / 1920.0, 155.0 / 1920.0), 40.0 / 1920.0, "494")) {
            std::cout << "wrong print" << '\n';
        }
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(392.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        DrawLine(input, glm::vec4(0.4, 0.4, 0.4, 0.7), glm::vec2(411.0 / 1920.0, 193.0 / 1920.0), glm::vec2(411.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(412.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(411.0 / 1920.0, 1764.0 / 1920.0), 40.0 / 1920.0, "29")) {
            std::cout << "wrong print" << '\n';
        }

        //weight
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(685.0 / 1920.0, 97.0 / 1920.0), 40.0 / 1920.0, "weight")) {
            std::cout << "wrong print" << '\n';
        }
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(685.0 / 1920.0, 155.0 / 1920.0), 40.0 / 1920.0, "5493")) {
            std::cout << "wrong print" << '\n';
        }
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(666.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        DrawLine(input, glm::vec4(0.4, 0.4, 0.4, 0.7), glm::vec2(685.0 / 1920.0, 193.0 / 1920.0), glm::vec2(685.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(686.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(685.0 / 1920.0, 1764.0 / 1920.0), 40.0 / 1920.0, "1260")) {
            std::cout << "wrong print" << '\n';
        }

        //horsepower
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(959.0 / 1920.0, 97.0 / 1920.0), 40.0 / 1920.0, "horsepower")) {
            std::cout << "wrong print" << '\n';
        }
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(959.0 / 1920.0, 155.0 / 1920.0), 40.0 / 1920.0, "249")) {
            std::cout << "wrong print" << '\n';
        }
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(940.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        DrawLine(input, glm::vec4(0.4, 0.4, 0.4, 0.7), glm::vec2(959.0 / 1920.0, 193.0 / 1920.0), glm::vec2(959.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(960.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(959.0 / 1920.0, 1764.0 / 1920.0), 40.0 / 1920.0, "27")) {
            std::cout << "wrong print" << '\n';
        }

        //accelerator
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1233.0 / 1920.0, 97.0 / 1920.0), 40.0 / 1920.0, "acceleration(0-60mph)")) {
            std::cout << "wrong print" << '\n';
        }
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1233.0 / 1920.0, 155.0 / 1920.0), 40.0 / 1920.0, "27")) {
            std::cout << "wrong print" << '\n';
        }
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(1214.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        DrawLine(input, glm::vec4(0.4, 0.4, 0.4, 0.7), glm::vec2(1233.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1233.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(1234.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1233.0 / 1920.0, 1764.0 / 1920.0), 40.0 / 1920.0, "6")) {
            std::cout << "wrong print" << '\n';
        }

        //mileage
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1507.0 / 1920.0, 97.0 / 1920.0), 40.0 / 1920.0, "mileage")) {
            std::cout << "wrong print" << '\n';
        }
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1507.0 / 1920.0, 155.0 / 1920.0), 40.0 / 1920.0, "51")) {
            std::cout << "wrong print" << '\n';
        }
        DrawFilledRect(input, glm::vec4(0.7, 0.6, 0.6, 0.95), glm::vec2(1488.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        DrawLine(input, glm::vec4(0.4, 0.4, 0.4, 0.7), glm::vec2(1507.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1507.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawFilledRect(input, glm::vec4(0.7, 0.6, 0.6, 0.95), glm::vec2(1508.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1507.0 / 1920.0, 1764.0 / 1920.0), 40.0 / 1920.0, "5")) {
            std::cout << "wrong print" << '\n';
        }

        //year
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1781.0 / 1920.0, 97.0 / 1920.0), 40.0 / 1920.0, "year")) {
            std::cout << "wrong print" << '\n';
        }
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1781.0 / 1920.0, 155.0 / 1920.0), 40.0 / 1920.0, "84")) {
            std::cout << "wrong print" << '\n';
        }
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(1762.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        DrawLine(input, glm::vec4(0.4, 0.4, 0.4, 0.7), glm::vec2(1781.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1781.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawFilledRect(input, glm::vec4(0.5, 0.5, 0.5, 0.2), glm::vec2(1782.0 / 1920.0, 193.0 / 1920.0), glm::vec2(18.0 / 1920.0, 1530.0 / 1920.0));
        if (! PrintText(input, glm::vec4(0, 0, 0, 0.7), glm::vec2(1781.0 / 1920.0, 1764.0 / 1920.0), 40.0 / 1920.0, "68")) {
            std::cout << "wrong print" << '\n';
        }
        

        std::vector<Car> data1(data);
        sort(data1.begin(), data1.end(), cmp);
        int num = data.size();
        glm::vec4 red(0.8, 0.4, 0.4, 0.8);
        glm::vec4 blue(0.4, 0.4, 0.8, 0.8);
        for (int i = 0; i < num; i++) {
            double cy = data1[i].cylinders;
            double dis = data1[i].displacement;
            double wei = data1[i].weight;
            double horse = data1[i].horsepower;
            double acce  = data1[i].acceleration;
            double milea = data1[i].mileage;
            double year  = data1[i].year;
            float xx1    = 0.0;
            float yy1    = 0.0;
            float xx2    = 0.0;
            float yy2    = 0.0;
            float xx3    = 0.0;
            float yy3    = 0.0;
            float xx4    = 0.0;
            float yy4    = 0.0;
            float xx5    = 0.0;
            float yy5    = 0.0;
            float xx6    = 0.0;
            float yy6    = 0.0;
            float xx7    = 0.0;
            float yy7    = 0.0;
            xx1          = 137.0 / 1920.0;
            yy1          = ((9 - cy) * 1.0 / 7.0 * 1530.0 + 193.0) / 1920.0;
            xx2          = 411.0 / 1920.0;
            yy2          = ((494 - dis) * 1.0 / 465.0 * 1530.0 + 193.0) / 1920.0;
            xx3          = 685.0 / 1920.0;
            yy3          = ((5493 - wei) * 1.0 / 4233.0 * 1530.0 + 193.0) / 1920.0;
            xx4          = 959.0 / 1920.0;
            yy4          = ((249 - horse) * 1.0 / 222.0 * 1530.0 + 193.0) / 1920.0;
            xx5          = 1233.0 / 1920.0;
            yy5          = ((27 - acce) * 1.0 / 21.0 * 1530.0 + 193.0) / 1920.0;
            xx6          = 1507.0 / 1920.0;
            yy6          = ((51 - milea) * 1.0 / 46.0 * 1530.0 + 193.0) / 1920.0;
            xx7          = 1781.0 / 1920.0;
            yy7          = ((84 - year) * 1.0 / 16.0 * 1530.0 + 193.0) / 1920.0;
            glm::vec4 col((red.r + i * 2.0 / num * 1.0 * (blue.r - red.r)), (red.g + i * 1.0 / num * 1.0 * (blue.g - red.g)), (red.b + i * 1.0 / num * 1.0 * (blue.b - red.b)), 0.8);
            
            DrawLine(input, col, glm::vec2(xx1, yy1), glm::vec2(xx2, yy2), 1.0);
            DrawLine(input, col, glm::vec2(xx2, yy2), glm::vec2(xx3, yy3), 1.0);
            DrawLine(input, col, glm::vec2(xx3, yy3), glm::vec2(xx4, yy4), 1.0);
            DrawLine(input, col, glm::vec2(xx4, yy4), glm::vec2(xx5, yy5), 1.0);
            DrawLine(input, col, glm::vec2(xx5, yy5), glm::vec2(xx6, yy6), 1.0);
            DrawLine(input, col, glm::vec2(xx6, yy6), glm::vec2(xx7, yy7), 1.0);
            
        }
        
        DrawLine(input, glm::vec4(1.0), glm::vec2(117.0 / 1920.0, 193.0 / 1920.0), glm::vec2(117.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(157.0 / 1920.0, 193.0 / 1920.0), glm::vec2(157.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(391.0 / 1920.0, 193.0 / 1920.0), glm::vec2(391.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(431.0 / 1920.0, 193.0 / 1920.0), glm::vec2(431.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(665.0 / 1920.0, 193.0 / 1920.0), glm::vec2(665.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(705.0 / 1920.0, 193.0 / 1920.0), glm::vec2(705.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(939.0 / 1920.0, 193.0 / 1920.0), glm::vec2(939.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(979.0 / 1920.0, 193.0 / 1920.0), glm::vec2(979.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(1213.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1213.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(1253.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1253.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(1487.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1487.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(1527.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1527.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(1761.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1761.0 / 1920.0, 1722.0 / 1920.0), 2.0);
        DrawLine(input, glm::vec4(1.0), glm::vec2(1801.0 / 1920.0, 193.0 / 1920.0), glm::vec2(1801.0 / 1920.0, 1722.0 / 1920.0), 2.0);

        return true;
    }

    void LIC(ImageRGB & output, Common::ImageRGB const & noise, VectorField2D const & field, int const & step) {
        // your code here
    }
}; // namespace VCX::Labs::Visualization