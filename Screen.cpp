#include "Screen.h"

vec2 ToScreenCoords(vec2 coords, int width, int height) {
    return vec2{
        ((coords.x + 1) / 2) * width,
        ((coords.y + 1) / 2) * height
    };
}