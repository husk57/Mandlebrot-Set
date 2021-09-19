//
//  main.cpp
//  Mandelbrot
//
//  Created by joshua on 9/18/21.
//

#include <iostream>
#include <cmath>
#include "tgaimage.h"

int width = 1024, height = 1024;

void createImg(int nmbr, double timeInterval) {
    TGAImage image(width, height, TGAImage::RGB);
    for (int y=0; y<width; y++) {
        for (int x=0; x<height; x++) {
            long double col[3] = {
                           0., 0., 0.
                           };
            long  double s[3] = {0.0, 0.0, 0.0};
            long double e[3] = {0.0, 1.0, 0.5};
            long double c[2] = {(long double)x/width, (long double)y/height};
            c[0] -= 0.5;
            c[1] -= 0.5;
            c[0] *= 3.5;
            c[1] *= 2.5;
            c[0] += -0.75;
            
            long double p[2] = {0.50011701007999932,0.54011501011999605};
            p[0] -= 0.5;
            p[1] -= 0.5;
            p[0] *= 3.5;
            p[1] *= 2.5;
            p[0] += -0.75;
            
            long double scale = pow(2.718, -timeInterval);
            long double depth = 0.0000000000000005;
            scale = depth + (1.-depth) * scale;
            
            c[0] = c[0] - p[0];
            c[1] = c[1] - p[1];
            c[0] *= depth;
            c[1] *= depth;
            c[0] = c[0] + p[0];
            c[1] = c[1] + p[1];
            
            long double z[3];
            std::copy(std::begin(c), std::end(c), std::begin(z));
            
            for (int i = 0; i<700; i++) {
                long double d = z[1];
                long double temp = d;
                d = z[0]*d + z[0]*d;
                long double outTemp[2] = {z[0]*z[0]+(temp*temp*-1.), d};
                outTemp[0] += c[0];
                outTemp[1] += c[1];
                std::copy(std::begin(outTemp), std::end(outTemp), std::begin(z));
                long double length = sqrt((z[0]*z[0])+(z[1]*z[1]));
                  if (length > 2.) {
                      col[0] = s[0] + (e[0]-s[0]) * (long double)(i)/50.;
                      col[1] = s[1] + (e[1]-s[1]) * (long double)(i)/50.;
                      col[2] = s[2] + (e[2]-s[2]) * (long double)(i)/50.;
                      break;
                  }
            }
            
            if (col[0] < 0.001) {
                col[0] = 0.;
            }
            if (col[1] < 0.001) {
                col[1] = 0.;
            }
            if (col[2] < 0.001) {
                col[2] = 0.;
            }
            image.set(x,y, TGAColor((int)(col[0] * 255.), (int)(col[1] * 255.), (int)(col[2] * 255.)));
        }
    }
    image.write_tga_file(std::to_string(nmbr) + ".tga");
}

int main(int argc, const char * argv[]) {
    long double time = 0.0;
    for (int i=0; i<1; i++)  {
        time += 1./2.;
        createImg(i, time);
        std::cout << i << std::endl;
    }
    std::cout << "done\n";
    return 0;
}
/*
 x0 := scaled x coordinate of pixel (scaled to lie in the Mandelbrot X scale (-2.5, 1))
     y0 := scaled y coordinate of pixel (scaled to lie in the Mandelbrot Y scale (-1, 1))
     x := 0.0
     y := 0.0
     iteration := 0
     max_iteration := 1000
     // Here N = 2^8 is chosen as a reasonable bailout radius.

     while x*x + y*y ≤ (1 << 16) and iteration < max_iteration do
         xtemp := x*x - y*y + x0
         y := 2*x*y + y0
         x := xtemp
         iteration := iteration + 1

     // Used to avoid floating point issues with points inside the set.
     if iteration < max_iteration then
         // sqrt of inner term removed using log simplification rules.
         log_zn := log(x*x + y*y) / 2
         nu := log(log_zn / log(2)) / log(2)
         // Rearranging the potential function.
         // Dividing log_zn by log(2) instead of log(N = 1<<8)
         // because we want the entire palette to range from the
         // center to radius 2, NOT our bailout radius.
         iteration := iteration + 1 - nu

     color1 := palette[floor(iteration)]
     color2 := palette[floor(iteration) + 1]
     // iteration % 1 = fractional part of iteration.
     color := linear_interpolate(color1, color2, iteration % 1)
     plot(Px, Py, color)
 */
/*
#version 410 core

uniform float time;
uniform vec2 resolution;
vec2 uv = gl_FragCoord.xy/resolution;
out vec4 FragColor;

void main()
{
    dvec2 c = uv;
    dvec3 s = vec3(0.0, 0.0, 0.0);
    dvec3 e = vec3(0.0, 1.0, 0.5);
    c -= vec2(0.5, 0.5);
    c *= vec2(3.5, 2.5);
    c += vec2(-.75, 0.);
    dvec2 p = vec2(0.500117,0.540115);
    p -= vec2(0.5, 0.5);
    p *= vec2(3.5, 2.5);
    p += vec2(-.75, 0.);
    double scale = pow(2.718, -time);
    scale = 0.000001 + (1.-0.000001) * scale;
    c = c - p;
    c *= 1.;
    c = c + p;
    dvec3 col = s;
    dvec2 z = c;
    for (int i =0; i<500; i++) {
      double d = z.y;
      double temp = d;
      d = z.x*d + z.x*d;
      z = dvec2(z.x*z.x+(temp*temp*-1.), d)+c;
        if (length(z) > 2.) {
          col = s + (e-s) * double(i)/50.;
            break;
        }
    }
    FragColor = vec4(col,1.0);
}
*/
