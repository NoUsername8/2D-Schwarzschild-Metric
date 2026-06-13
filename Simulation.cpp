#include "math/Matrix.hpp"
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
#define SCALE 25
#define THICKNESS 0.98
#define boundary 15
#define TOLERANCE 1e-10

bool text = false;
bool tp = false;
bool cordGrid = false;

#include "Metric.h"
#include "ScreenHelper.h"
#include "math/Polar.hpp"
#include "math/Vector.hpp"

Vector<8> y[] = {
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0025}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.005}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0075}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.01}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0125}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.015}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0175}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.02}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0225}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.025}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0275}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.03}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0325}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.035}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0375}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.04}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0425}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.045}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.0475}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, 0.05}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0025}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.005}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0075}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.01}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0125}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.015}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0175}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.02}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0225}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.025}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0275}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.03}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0325}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.035}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0375}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.04}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0425}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.045}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.0475}),
    Vector<8>({10, -1, 0, 0, M_PI_2, 0, 0, -0.05}),
    Vector<8>({1.5 * rs, 0, 0, 0, M_PI_2, 0, M_PI_2, 1}),
    Vector<8>({1.25 * rs, 0, 0, 0, M_PI_2, 0, M_PI_2, 1}),
    Vector<8>({1.05 * rs, 0, 0, 0, M_PI_2, 0, M_PI_2, 1}),
};
constexpr int pointSize = sizeof(y) / sizeof(y[0]);

constexpr double timeSpan = 100;
constexpr double timeStep = 0.005;

constexpr int movementRate = 2;

Vector<3> variant1(Vector<3> y) { return y + timeStep * getYPrime(y); }

Vector<3> variant2(Vector<3> y) {

  double vals[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Matrix<3, 3> id3x3(vals);

  Vector<3> yNew = y;

  for (int i = 0; i < 200; i++) {

    // vec6 phi = (yNew - y) * (1 / timeStep) - getYPrime(yNew);
    // matrix6x6 phi_ = (1 / timeStep) * id6x6() - getYPrimePrime(yNew);
    Vector<3> phi = yNew - y - timeStep * getYPrime((yNew + y) * 0.5);
    Matrix<3, 3> phi_ =
        id3x3 - timeStep * 0.5 * getYPrimePrime((yNew + y) * 0.5);
    printf("%d PHI: %s \n", i, phi.toString().c_str());
    y = yNew;
    yNew = y - inv(phi_) * phi;
  }
  y = yNew;

  return y;
}

int main() {

  // Setup SDL
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    printf("Couldn't initialize SDL: %s\n", SDL_GetError());
    return EXIT_FAILURE;
  }
  SDL_Window *window =
      SDL_CreateWindow("Example: 0", SDL_WINDOWPOS_UNDEFINED,
                       SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, 0);
  if (!window) {
    printf("Failed to open %d x %d window: %s\n", SCREEN_WIDTH, SCREEN_HEIGHT,
           SDL_GetError());
    return EXIT_FAILURE;
  }
  SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
  SDL_RenderClear(renderer);
  SDL_Event event;
  bool up = false;
  bool down = false;
  bool left = false;
  bool right = false;
  bool space = false;
  bool paused = false;
  bool sPressed = false;
  Vector<2> windowPos({0, 0});

  for (int i = 0; i < pointSize; i++) {
    // Initial Normalisation
    double velocityMagnitude =
        std::sqrt(y[i](1) * y[i](1) + y[i](3) * y[i](3) + y[i](5) * y[i](5) +
                  y[i](7) * y[i](7)) != 0
            ? std::sqrt(y[i](1) * y[i](1) + y[i](3) * y[i](3) +
                        y[i](5) * y[i](5) + y[i](7) * y[i](7))
            : 1;
    y[i](1) = (y[i](1) / velocityMagnitude);
    y[i](3) = (y[i](3) / velocityMagnitude);
    y[i](5) = (y[i](5) / velocityMagnitude);
    y[i](7) = (y[i](7) / velocityMagnitude);
  }

  std::list<Vector<2>> trajectory[pointSize];
  bool sim[pointSize];
  for (int i = 0; i < pointSize; i++) {
    sim[i] = true;
  }

  for (int t = 0; t <= timeSpan; t += timeStep) {
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);

    // Event Handler
    do {
      while (SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_QUIT:
          t = timeSpan + 1;
          break;
        case SDL_KEYDOWN:
          switch (event.key.keysym.scancode) {
          case SDL_SCANCODE_UP:
            up = true;
            break;
          case SDL_SCANCODE_DOWN:
            down = true;
            break;
          case SDL_SCANCODE_LEFT:
            left = true;
            break;
          case SDL_SCANCODE_RIGHT:
            right = true;
            break;
          case SDL_SCANCODE_SPACE:
            space = true;
            break;
          }
          break;
        case SDL_KEYUP:
          switch (event.key.keysym.scancode) {
          case SDL_SCANCODE_UP:
            up = false;
            break;
          case SDL_SCANCODE_DOWN:
            down = false;
            break;
          case SDL_SCANCODE_LEFT:
            left = false;
            break;
          case SDL_SCANCODE_RIGHT:
            right = false;
            break;
          case SDL_SCANCODE_SPACE:
            space = false;
            sPressed = false;
            break;
          }
          break;
        default:
          break;
        }
        if (space && !sPressed) {
          sPressed = true;
          paused = !paused;
        }
        if (up) {
          windowPos(1) += movementRate;
        }
        if (down) {
          windowPos(1) -= movementRate;
        }
        if (left) {
          windowPos(0) += movementRate;
        }
        if (right) {
          windowPos(0) -= movementRate;
        }
      }
    } while (paused);

    // Draw Coordinate Grid
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      for (int y = 0; y < SCREEN_HEIGHT; y++) {
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // WHITE
        Vector<2> pos = TransformToSimulationCoords(x, y, windowPos);
        Vector<2> polarPosition = cartesianToPolar(Vector<2>({pos(0), pos(1)}));
        double r = polarPosition(0);
        double theta = polarPosition(1);
        if (std::fmod(r, 1) >= THICKNESS && cordGrid) {
          SDL_RenderDrawPoint(renderer, x, y);
        }
        if (r >= rs && r <= rs + THICKNESS / SCALE) {
          SDL_SetRenderDrawColor(renderer, 255, 0, 255,
                                 255); // PURPLE
          SDL_RenderDrawPoint(renderer, x, y);
        }
      }
    }

    if (cordGrid) {
      SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // WHITE
      Vector<2> originOnScreen = TransformToScreenCoords(0, 0, windowPos);
      for (int i = 0; i < 8; i++) {
        int size = 500000;
        Vector<2> screenEndPoint =
            TransformToScreenCoords(size * std::cos(i * M_PI_4),
                                    size * std::sin(i * M_PI_4), windowPos);
        SDL_RenderDrawLine(renderer, originOnScreen(0), originOnScreen(1),
                           screenEndPoint(0), screenEndPoint(1));
      }
    }
    for (int i = 0; i < pointSize; i++) {

      // Velocity
      double vT = y[i](3);
      double vR = y[i](1);
      double vA = y[i](7);
      double vP = y[i](5);

      double pointR = y[i](0);
      double pointT = y[i](6);

      double velocity = std::sqrt(vT * vT + vR * vR + vA * vA + vP * vP);

      if (pointR < rs + TOLERANCE) {

        sim[i] = false;
      }

      // Display Basis
      /*Basis basis = getPolarBasis(pointR, pointT);
      Vector<2> probepoint{pointR * -std::cos(pointT), pointR *
      std::sin(pointT)}; Vector<2> rEndpoint{basis.e1.x + probepoint.x,
      basis.e1.y + probepoint.y}; Vector<2> tEndpoint{basis.e2.x +
      probepoint.x, basis.e2.y + probepoint.y}; Vector<2> probeOnScreen =
          TransformToScreenCoords(probepoint.x, probepoint.y, windowPos);
      Vector<2> rPointOnScreen =
          TransformToScreenCoords(rEndpoint.x, rEndpoint.y, windowPos);
      Vector<2> tPointOnScreen =
          TransformToScreenCoords(tEndpoint.x, tEndpoint.y, windowPos);
      SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); // RED
      SDL_RenderDrawLine(renderer, probeOnScreen.x, probeOnScreen.y,
                         rPointOnScreen.x, rPointOnScreen.y);
      SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255); // GREEN
      SDL_RenderDrawLine(renderer, probeOnScreen.x, probeOnScreen.y,
                         tPointOnScreen.x, tPointOnScreen.y);*/

      if (sim[i]) {
        Vector<3> y3 = Vector<3>({y[i](0), y[i](1), y[i](7)});
        y3 = variant2(y3);
        y[i](0) = y3(0);
        y[i](1) = y3(1);
        y[i](6) += timeStep * y3(2);
        y[i](7) = y3(2);

        // Add point to Trajectory
        Vector<2> cartesianPosition =
            polarToCartesian(Vector<2>({pointR, pointT}));
        trajectory[i].push_back(cartesianPosition);
      }
      // Draw Trajectory
      SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255); // BLUE
      for (Vector<2> point : trajectory[i]) {
        Vector<2> screenPoint =
            TransformToScreenCoords(point(0), point(1), windowPos);
        SDL_RenderDrawPoint(renderer, screenPoint(0), screenPoint(1));
      }

      // Iterative Normalisation
      double velocityMagnitude =
          std::sqrt(y[i](1) * y[i](1) + y[i](3) * y[i](3) + y[i](5) * y[i](5) +
                    y[i](7) * y[i](7)) != 0
              ? std::sqrt(y[i](1) * y[i](1) + y[i](3) * y[i](3) +
                          y[i](5) * y[i](5) + y[i](7) * y[i](7))
              : 1;
      y[i](1) = (y[i](1) / velocityMagnitude);
      y[i](3) = (y[i](3) / velocityMagnitude);
      y[i](5) = (y[i](5) / velocityMagnitude);
      y[i](7) = (y[i](7) / velocityMagnitude);

      // Teleport
      if (tp) {
        Vector<2> cartPos = polarToCartesian(Vector<2>({y[i](0), y[i](6)}));
        if (y[i](0) <= 1 && false) {
          t = timeSpan + 1;
        }
        if (y[i](0) >= boundary) {
          y[i](0) = boundary - 1;
          y[i](6) = cartesianToPolar(Vector<2>({-cartPos(1), -cartPos(2)}))(1);
          y[i](1) = -y[i](1);
          y[i](7) = -y[i](7);
        }
      }

      // Text
      if (text && i == 0) {
        TTF_Init();
        TTF_Font *Sans =
            TTF_OpenFont("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 24);
        SDL_Color Gray = {255, 255, 255};
        std::string out = "V-Time: " + std::to_string(y[i](3)) +
                          " V-Radial: " + std::to_string(y[i](1)) +
                          " V-Azimuth: " + std::to_string(y[i](7)) +
                          " V-Polar: " + std::to_string(y[i](5)) +
                          " V-Total: " + std::to_string(velocity);
        SDL_Surface *surfaceMessage =
            TTF_RenderText_Solid(Sans, out.c_str(), Gray);
        SDL_Texture *Message =
            SDL_CreateTextureFromSurface(renderer, surfaceMessage);
        SDL_Rect Message_rect;
        Message_rect.w = 1000;
        Message_rect.h = 2 * 24;
        Message_rect.x = SCREEN_WIDTH - Message_rect.w;
        Message_rect.y = SCREEN_HEIGHT - Message_rect.h;
        SDL_RenderCopy(renderer, Message, NULL, &Message_rect);
        SDL_FreeSurface(surfaceMessage);
        SDL_DestroyTexture(Message);
        TTF_CloseFont(Sans);
        TTF_Quit();
      }
    }

    // SDL stuff
    SDL_RenderPresent(renderer);
    SDL_Delay(1000 * timeStep);
  }

  // Do SDL stuff
  SDL_RenderPresent(renderer);
  SDL_DestroyWindow(window);
  SDL_DestroyRenderer(renderer);
  SDL_Quit();
  return EXIT_SUCCESS;
}
