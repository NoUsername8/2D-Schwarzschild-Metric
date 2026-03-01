#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <list>

#include "GeometryHelper.h"
#include "Screen.h"

#define _USE_MATH_DEFINES
#include <cmath>

#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
#define SCALE 120

#define LIMIT 10000
#define STEP 0.1
constexpr double GRIDSIZE = 1;

struct Probe{
  vec2 position;
  double distance;
  int movingBasis;
  Basis basis;
  bool first;
};

int main() {

  //Setup SDL
  if (SDL_Init(SDL_INIT_VIDEO) < 0){
    printf("Couldn't initialize SDL: %s\n", SDL_GetError());
    return EXIT_FAILURE;
  }
  SDL_Window *window = SDL_CreateWindow("Example: 0", SDL_WINDOWPOS_UNDEFINED,
					SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, 0);
  if (!window){
    printf("Failed to open %d x %d window: %s\n", SCREEN_WIDTH, SCREEN_HEIGHT, SDL_GetError());
    return EXIT_FAILURE;
  }
  SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
  SDL_RenderClear(renderer);
  SDL_Event event;

  //Setup Simulation
  bool running = true;
  int index = 0;
  std::vector<Probe> probes;
  std::list<vec2> greenPoints;
  std::list<vec2> redPoints;
  int pSize = 4;
  probes.push_back(Probe{vec2{0, 0}, 0, 1, true});
  probes.push_back(Probe{vec2{0, 0}, 0, -1, true});
  probes.push_back(Probe{vec2{0, 0}, 0, 2, true});
  probes.push_back(Probe{vec2{0, 0}, 0, -2, true});

  //Timestep
  for(int i = 0; i < LIMIT; i++) {

    //For every probe
    for(Probe p : probes) {

      //Spawn new probe
      if(fmod(p.distance, GRIDSIZE) == 0 && p.first) {
        int direction = std::abs(p.movingBasis) == 1 ? 2 : 1;
        probes.push_back(Probe{p.position, 0, direction, false});
        probes.push_back(Probe{p.position, 0, -direction, false});
      }

      //Move Probe
      switch(p.movingBasis) {
        case 1:
          p.position = p.position + p.basis.e1 * STEP;
          break;
        case -1:
          p.position = p.position + p.basis.e1 * -STEP;
          break;
        case 2:
          p.position = p.position + p.basis.e2 * STEP;
          break;
        case -2:
          p.position = p.position + p.basis.e2 * -STEP;
          break;
      }

      //Update Probe
      p.basis = getBasis(p.position);
      p.distance += STEP;

      //Add Point to List
      if(std::abs(p.movingBasis) == 1) {
        redPoints.push_back(p.position);
      } else {
        greenPoints.push_back(p.position);
      }
    }

    //Render
    SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
    for(vec2 v : greenPoints) {
      vec2 screen = ToScreenCoords(v, SCREEN_WIDTH, SCREEN_HEIGHT);
      SDL_RenderDrawPoint(renderer, screen.x, screen.y);
    }
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    for(vec2 v : redPoints) {
      vec2 screen = ToScreenCoords(v, SCREEN_WIDTH, SCREEN_HEIGHT);
      SDL_RenderDrawPoint(renderer, screen.x, screen.y);
    }
    SDL_RenderPresent(renderer);

    //Update Array Size
    pSize = sizeof(probes) / sizeof(Probe);
  }

  //Wait
  std::cout << "DONE!" << std::endl;
  while(true);

  //Quit
  SDL_DestroyWindow(window);
  SDL_DestroyRenderer(renderer);
  SDL_Quit();
  return EXIT_SUCCESS;
}

