/* FT-GPI 
 * Copyright (C) 2022 Rodrigo Canovas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/ .
 * */


#include <vector>
#include <cmath>

#ifndef FDGPI_TABLES_H
#define FDGPI_TABLES_H

namespace fdgpi {

    std::string rowcode = "ACDEFGHIKLMNPQRSTVWY"; //for the matrices

    const int n_col = 21;
    const int n_row = 20;

    /*Center of TM-helix (inside to outside)
     * matrix of size 20x21   (20 -> one for each amino acid)
     * columns -> position from TMS-center
     *  rows  -> amino acid order as ACDEFGHIKLMNPQRSTVWY
     */
    const std::vector<std::vector<double> > io_center = {
            {0.00, 0.00,  0.00,   0.00,  0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,  0.00, 0.00,   0.00,   0.00,  0.00,  0.00,   0.00,  0.00,  0.00,  0.00,  0.00},
            {0.00, 56.22, 63.83,  56.55, 56.78, 50.81, 62.05,  61.56,  53.03,  50.50,  90.20, 75.05, 66.69,  59.31,  58.52, 64.03, 54.84,  70.12, 60.23, 64.98, 51.71, 43.31},
            {0.00, 7.63,  12.32,  14.62, 17.69, 13.12, 13.58,  15.09,  12.44,  12.18,  19.59, 26.94, 26.83,  22.94,  16.19, 13.03, 11.59,  14.18, 18.36, 8.41,  11.71, 7.33},
            {0.00, 1.50,  5.53,   3.13,  5.12,  1.65,  9.22,   9.16,   9.97,   7.94,   6.60,  3.97,  1.65,   3.64,   2.83,  6.50,  0.21,   6.90,  2.89,  5.04,  3.81,  4.00},
            {0.00, 1.67,  2.97,   3.50,  3.10,  2.11,  3.75,   2.12,   4.00,   3.50,   0.11,  1.54,  3.33,   3.45,   8.55,  1.98,  10.03,  7.66,  4.02,  7.36,  1.09,  2.05},
            {0.00, 36.67, 45.11,  62.78, 45.64, 51.14, 53.22,  48.46,  65.50,  66.48,  60.30, 58.78, 51.14,  55.38,  46.57, 53.64, 47.74,  58.07, 42.99, 44.93, 58.14, 48.39},
            {0.00, 33.83, 35.91,  38.16, 36.76, 46.55, 50.70,  40.60,  41.82,  38.26,  38.43, 40.66, 40.50,  54.45,  39.37, 51.94, 56.37,  42.83, 40.84, 44.38, 51.44, 49.29},
            {0.00, 5.09,  6.30,   4.16,  5.78,  5.80,  3.42,   5.58,   4.21,   1.00,   3.67,  2.01,  4.83,   0.75,   3.25,  5.05,  10.50,  9.88,  3.75,  5.58,  5.75,  2.00},
            {0.00, 45.58, 67.14,  64.72, 79.95, 71.41, 65.55,  69.56,  72.82,  70.75,  74.20, 69.14, 71.26,  58.16,  73.63, 66.33, 59.31,  71.02, 76.40, 57.96, 55.24, 48.86},
            {0.00, 5.79,  1.70,   8.30,  4.62,  5.07,  1.00,   1.00,   3.00,   2.00,   3.00,  1.25,  1.00,   4.50,   2.00,  1.00,  6.00,   1.92,  4.87,  7.00,  1.50,  1.53},
            {0.00, 73.51, 104.45, 99.52, 88.65, 93.90, 114.73, 106.01, 117.54, 106.47, 85.43, 94.79, 122.26, 101.94, 91.46, 96.13, 107.92, 85.98, 73.36, 91.43, 92.36, 76.41},
            {0.00, 16.54, 22.80,  24.90, 33.16, 22.61, 15.64,  19.39,  13.82,  20.57,  19.12, 23.47, 22.39,  20.70,  24.21, 21.49, 18.96,  12.84, 31.23, 23.70, 29.09, 16.94},
            {0.00, 7.44,  9.91,   14.42, 14.06, 15.48, 11.42,  6.27,   9.73,   9.76,   11.32, 7.11,  3.53,   8.18,   12.82, 9.72,  13.23,  8.53,  14.07, 9.43,  9.72,  10.11},
            {0.00, 8.99,  14.01,  8.19,  8.04,  9.71,  13.76,  13.38,  7.75,   10.84,  8.79,  11.91, 7.54,   11.75,  18.90, 15.67, 22.76,  33.07, 30.79, 24.90, 13.90, 13.27},
            {0.00, 3.11,  7.51,   8.10,  9.29,  2.99,  10.34,  9.83,   7.65,   10.80,  7.07,  10.83, 8.22,   6.93,   6.37,  8.10,  7.01,   10.84, 6.50,  8.52,  9.02,  5.19},
            {0.00, 1.20,  3.03,   5.33,  5.19,  4.08,  4.24,   2.50,   1.25,   1.36,   1.50,  0.50,  4.17,   1.75,   3.00,  3.00,  2.50,   3.36,  3.00,  8.09,  5.68,  1.06},
            {0.00, 14.77, 30.12,  28.53, 37.06, 45.56, 24.71,  39.87,  31.55,  33.63,  46.58, 40.50, 35.65,  38.41,  50.45, 35.50, 26.92,  30.12, 34.85, 33.97, 28.44, 26.55},
            {0.00, 12.00, 23.80,  35.63, 26.46, 37.10, 27.14,  39.74,  29.40,  26.06,  25.38, 28.21, 21.55,  39.66,  31.72, 40.34, 33.98,  33.29, 24.81, 31.47, 23.54, 25.38},
            {0.00, 56.55, 63.88,  68.47, 75.25, 67.42, 70.65,  65.20,  62.34,  79.23,  66.26, 65.09, 63.71,  68.46,  70.35, 56.01, 58.67,  53.80, 65.92, 62.42, 54.49, 48.86},
            {0.00, 5.75,  10.55,  10.31, 14.06, 16.44, 13.88,  19.43,  27.22,  17.57,  12.88, 17.29, 16.00,  14.20,  13.84, 15.87, 10.25,  7.58,  17.21, 13.06, 15.26, 15.07},
            {0.00, 19.51, 15.99,  16.51, 12.37, 18.25, 14.20,  11.46,  12.18,  18.30,  6.77,  8.18,  14.97,  12.65,  13.18, 21.88, 27.43,  19.22, 22.92, 26.08, 36.92, 39.03}
    };

    /*Center of TM-helix (outside to inside)
    * matrix of size 20x21   (20 -> one for each amino acid)
    * columns -> position from TMS-center
    *  rows  -> amino acid order as ACDEFGHIKLMNPQRSTVWY
    */
    const std::vector<std::vector<double> > oi_center = {
            {0.00, 0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00},
            {0.00, 64.54, 58.99,  72.81,  79.11,  72.30,  79.80,  87.32,  87.37,  87.75,  81.26,  69.69,  95.78,  74.08,  75.32,  61.43,  66.79,  66.99,  62.01,  70.68,  60.23,  57.11},
            {0.00, 2.79,  4.62,   10.11,  15.80,  19.29,  17.48,  13.95,  25.93,  16.50,  18.01,  21.73,  19.31,  17.77,  21.84,  20.20,  22.83,  14.94,  19.52,  17.50,  20.75,  23.31},
            {0.00, 4.63,  5.10,   2.14,   6.99,   7.16,   8.74,   4.07,   5.05,   6.25,   6.03,   2.50,   2.61,   2.00,   2.25,   4.88,   2.50,   2.00,   3.67,   4.00,   3.06,   2.72},
            {0.00, 5.38,  5.87,   3.92,   5.36,   12.37,  7.48,   6.42,   2.20,   3.25,   5.26,   4.29,   6.52,   8.69,   9.96,   5.82,   4.41,   3.17,   12.06,  8.03,   8.16,   4.56},
            {0.00, 49.46, 48.03,  61.90,  46.07,  65.27,  61.62,  59.28,  71.32,  74.85,  63.80,  67.43,  56.76,  59.54,  51.83,  63.49,  62.80,  51.98,  57.24,  63.33,  71.92,  72.01},
            {0.00, 37.47, 48.20,  56.71,  63.42,  57.02,  68.19,  66.20,  74.68,  73.53,  62.81,  56.85,  49.80,  60.49,  60.43,  64.78,  56.70,  49.52,  55.92,  55.75,  31.98,  31.15},
            {0.00, 4.58,  9.58,   3.14,   7.00,   5.58,   3.75,   1.13,   4.83,   5.00,   2.51,   6.70,   1.90,   2.73,   1.50,   3.59,   1.59,   5.71,   1.58,   9.72,   4.53,   5.09},
            {0.00, 76.94, 81.49,  93.25,  95.98,  85.71,  86.77,  94.03,  82.35,  87.49,  97.25,  93.86,  89.50,  88.49,  91.90,  82.56,  86.35,  94.59,  101.14, 65.92,  84.86,  67.32},
            {0.00, 3.13,  3.96,   3.34,   2.48,   5.28,   4.79,   2.87,   3.85,   3.20,   2.71,   3.38,   2.39,   5.36,   3.56,   1.60,   6.59,   9.83,   7.41,   8.15,   9.20,   1.06},
            {0.00, 86.53, 107.26, 131.25, 128.92, 110.43, 112.79, 118.47, 121.83, 132.11, 134.37, 132.62, 138.71, 139.10, 149.01, 162.76, 142.80, 146.29, 131.53, 122.70, 124.00, 112.73},
            {0.00, 19.37, 34.04,  23.04,  18.65,  25.87,  21.86,  25.11,  25.82,  20.93,  24.23,  23.46,  25.86,  26.79,  23.96,  24.00,  22.18,  33.36,  25.93,  22.95,  17.15,  25.34},
            {0.00, 5.98,  8.44,   13.31,  8.70,   6.31,   8.71,   5.50,   10.71,  13.95,  16.33,  9.43,   11.56,  10.53,  20.32,  23.93,  21.68,  18.01,  11.03,  14.61,  9.78,   9.38},
            {0.00, 19.53, 26.97,  20.41,  18.42,  26.21,  17.05,  13.00,  17.61,  7.69,   15.87,  14.37,  28.66,  23.52,  26.74,  23.34,  28.56,  19.56,  10.86,  12.78,  9.77,   19.50},
            {0.00, 6.07,  7.42,   13.56,  12.00,  2.34,   11.18,  3.20,   3.00,   5.67,   4.83,   2.31,   4.73,   3.67,   4.33,   5.32,   5.37,   6.76,   7.44,   6.33,   9.40,   3.36},
            {0.00, 3.07,  6.08,   8.19,   6.16,   9.63,   8.10,   3.06,   3.27,   2.85,   4.56,   4.56,   1.32,   4.20,   2.56,   3.00,   4.64,   11.49,  7.20,   13.07,  5.74,   6.71},
            {0.00, 20.75, 37.03,  39.40,  48.76,  53.36,  56.34,  45.96,  43.70,  47.73,  42.91,  49.64,  55.59,  52.56,  42.97,  40.51,  39.83,  37.24,  39.84,  36.49,  42.33,  26.95},
            {0.00, 23.09, 28.63,  38.89,  35.50,  42.54,  32.60,  40.43,  39.06,  35.82,  38.05,  38.71,  35.16,  37.43,  35.02,  34.51,  36.51,  29.05,  38.53,  34.98,  29.11,  24.67},
            {0.00, 58.99, 98.02,  87.41,  84.01,  86.55,  95.13,  100.90, 90.45,  84.63,  94.42,  110.96, 92.77,  98.78,  90.64,  90.34,  81.65,  96.36,  91.62,  96.34,  86.32,  69.76},
            {0.00, 27.83, 35.01,  19.34,  20.58,  16.57,  11.26,  15.23,  6.65,   13.87,  10.38,  8.48,   7.45,   7.17,   8.56,   12.56,  14.83,  11.98,  14.35,  14.51,  18.48,  18.81},
            {0.00, 31.05, 30.86,  23.80,  29.23,  23.88,  22.56,  30.07,  17.51,  14.14,  11.62,  16.22,  10.82,  14.29,  14.50,  8.58,   27.60,  25.85,  34.42,  51.28,  53.52,  49.94}
    };

    /*N-terminus of TM-helix (inside to outside)
    * matrix of size 20x21   (20 -> one for each amino acid)
    * columns -> position from TMS-center
    *  rows  -> amino acid order as ACDEFGHIKLMNPQRSTVWY
    */
    const std::vector<std::vector<double> > io_nt = {
            {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
            {0.00, 35.25, 33.27, 33.54, 52.30, 26.50,  76.76,  63.82,  54.85,  63.23, 45.96,  64.82,  64.65,  57.72,  61.81,  83.74,  56.27,  68.55,  63.29, 65.44,  49.73, 69.97},
            {0.00, 6.83,  6.57,  12.94, 3.72,  7.09,   8.46,   17.66,  11.99,  17.32, 11.67,  16.65,  12.74,  10.65,  16.73,  14.59,  24.11,  30.03,  24.18, 16.31,  15.43, 14.25},
            {0.00, 26.28, 12.54, 12.63, 22.32, 24.53,  1.19,   4.00,   5.17,   4.29,  4.68,   5.98,   7.57,   5.58,   11.49,  9.65,   3.10,   1.00,   4.81,  3.67,   2.67,  2.51},
            {0.00, 17.30, 39.69, 13.54, 15.90, 29.97,  4.11,   6.00,   5.29,   2.17,  3.67,   4.40,   2.50,   1.94,   2.50,   1.50,   1.00,   3.58,   0.52,  7.81,   8.39,  11.46},
            {0.00, 29.45, 36.43, 28.90, 23.54, 20.21,  52.08,  62.29,  51.28,  47.28, 47.47,  57.80,  51.74,  60.52,  64.92,  58.14,  62.27,  47.53,  52.45, 54.87,  40.52, 63.86},
            {0.00, 45.83, 27.63, 39.87, 40.93, 22.84,  45.14,  32.45,  48.16,  48.45, 41.07,  34.33,  43.21,  39.85,  43.34,  37.53,  45.09,  43.86,  59.70, 44.43,  52.95, 50.38},
            {0.00, 17.79, 15.17, 17.27, 13.78, 24.26,  5.90,   8.05,   4.35,   5.66,  6.02,   6.50,   4.51,   3.38,   3.00,   2.67,   6.50,   0.25,   1.33,  3.05,   7.00,  9.47},
            {0.00, 29.64, 22.68, 45.13, 24.12, 11.78,  72.23,  70.95,  53.36,  78.85, 91.16,  62.54,  69.34,  73.38,  59.69,  66.46,  75.16,  76.48,  60.11, 64.12,  74.39, 56.40},
            {0.00, 54.67, 50.46, 43.82, 38.35, 76.72,  5.34,   3.49,   5.26,   8.07,  3.67,   1.00,   1.00,   3.00,   2.00,   1.00,   1.00,   5.50,   2.25,  2.00,   5.00,  4.25},
            {0.00, 45.56, 49.57, 73.54, 52.53, 28.63,  105.25, 110.59, 110.50, 82.26, 101.00, 103.93, 118.42, 106.31, 101.58, 107.29, 105.41, 113.91, 92.10, 109.79, 85.10, 84.34},
            {0.00, 20.26, 14.35, 10.93, 20.76, 11.72,  24.35,  24.39,  35.99,  24.33, 22.30,  22.90,  10.99,  18.25,  22.78,  17.73,  25.47,  18.05,  22.49, 25.17,  19.15, 13.39},
            {0.00, 24.52, 21.44, 28.08, 30.65, 25.74,  6.16,   11.25,  17.50,  11.53, 17.43,  9.26,   8.92,   12.12,  8.91,   10.76,  3.91,   3.83,   12.82, 8.27,   18.13, 8.96},
            {0.00, 21.28, 33.71, 28.14, 26.46, 18.01,  14.00,  9.44,   13.67,  8.54,  11.01,  8.72,   12.72,  12.50,  11.00,  11.84,  10.62,  2.41,   8.56,  15.39,  27.52, 20.52},
            {0.00, 24.87, 23.03, 12.10, 13.37, 22.29,  7.03,   6.85,   10.05,  8.20,  4.90,   11.08,  5.70,   12.45,  9.49,   10.09,  9.64,   6.07,   4.64,  4.99,   8.11,  10.62},
            {0.00, 65.36, 45.71, 51.35, 65.20, 114.55, 1.58,   3.76,   5.93,   4.53,  3.31,   1.11,   5.50,   4.67,   1.67,   0.17,   2.42,   1.50,   5.25,  4.50,   1.00,  1.11},
            {0.00, 40.23, 51.10, 35.44, 38.80, 51.31,  19.44,  28.81,  32.54,  39.26, 39.36,  32.00,  39.06,  38.08,  30.28,  35.61,  42.48,  51.94,  39.37, 38.67,  38.32, 44.34},
            {0.00, 23.71, 47.77, 26.25, 27.78, 30.48,  17.42,  29.66,  33.98,  35.41, 28.92,  45.14,  28.03,  27.06,  23.33,  36.62,  19.58,  24.24,  37.19, 34.29,  30.89, 32.19},
            {0.00, 35.56, 30.50, 42.30, 31.88, 25.40,  74.48,  64.55,  62.05,  65.29, 84.42,  80.65,  69.52,  71.06,  68.25,  61.98,  76.66,  57.42,  63.91, 66.14,  66.97, 50.08},
            {0.00, 10.54, 13.86, 14.33, 20.30, 3.67,   22.82,  14.61,  11.95,  17.28, 8.88,   11.00,  16.03,  19.02,  30.68,  15.72,  9.56,   18.04,  24.11, 12.08,  17.99, 11.19},
            {0.00, 15.16, 14.61, 19.98, 27.40, 14.40,  26.33,  17.47,  16.21,  18.13, 13.19,  10.28,  17.94,  12.53,  16.64,  6.98,   9.83,   15.87,  11.02, 8.10,   19.82, 26.80}
    };

    /*N-terminus of TM-helix (outside to inside)
    * matrix of size 20x21   (20 -> one for each amino acid)
    * columns -> position from TMS-center
    *  rows  -> amino acid order as ACDEFGHIKLMNPQRSTVWY
    */
    const std::vector<std::vector<double> > oi_nt = {
            {0.00, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,   0.00,   0.00,   0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00},
            {0.00, 59.90, 51.67, 71.84, 46.91, 30.97, 71.37,  67.31,  72.18,  67.38, 85.32,  76.99,  83.93,  92.89,  85.76,  87.06,  79.01,  88.77,  78.53,  84.38,  62.91,  65.97},
            {0.00, 20.54, 9.62,  2.64,  15.07, 6.69,  3.33,   7.37,   10.41,  12.41, 17.68,  17.66,  16.94,  19.38,  22.36,  20.83,  15.94,  25.69,  15.95,  19.51,  17.91,  22.33},
            {0.00, 29.04, 32.96, 29.27, 29.28, 54.08, 7.82,   1.87,   11.45,  3.49,  7.36,   8.42,   4.66,   2.07,   6.53,   3.05,   2.21,   5.65,   3.63,   1.50,   3.50,   3.00},
            {0.00, 42.88, 37.55, 34.49, 35.12, 54.94, 7.45,   5.09,   5.70,   7.82,  10.40,  9.85,   4.21,   2.74,   3.30,   3.40,   1.67,   8.37,   6.05,   5.62,   7.54,   11.55},
            {0.00, 36.98, 42.99, 31.09, 30.16, 9.78,  65.59,  59.36,  74.72,  60.29, 53.20,  60.76,  52.82,  63.31,  65.77,  69.93,  68.68,  67.62,  60.58,  49.50,  64.49,  62.14},
            {0.00, 58.54, 53.94, 60.03, 47.22, 42.39, 52.37,  60.41,  48.16,  68.29, 51.00,  65.60,  77.12,  67.85,  72.66,  58.34,  57.89,  57.79,  62.71,  63.04,  63.92,  48.58},
            {0.00, 13.26, 17.91, 23.54, 23.29, 29.56, 6.21,   9.40,   5.71,   9.62,  5.27,   3.92,   3.43,   7.17,   3.83,   4.50,   3.50,   4.80,   1.50,   4.07,   0.00,   5.40},
            {0.00, 30.96, 45.59, 38.79, 29.21, 21.28, 97.55,  74.14,  101.28, 93.19, 100.12, 94.18,  85.83,  81.68,  85.77,  94.62,  104.99, 96.24,  81.46,  97.28,  82.10,  85.56},
            {0.00, 38.70, 20.95, 38.67, 37.28, 74.78, 3.21,   4.40,   7.01,   4.62,  5.14,   1.53,   5.53,   1.80,   2.50,   4.53,   5.91,   1.94,   2.02,   6.96,   3.40,   5.52},
            {0.00, 61.84, 63.02, 68.70, 78.66, 40.99, 132.30, 135.72, 124.33, 91.97, 124.04, 131.18, 114.95, 112.00, 138.28, 133.09, 136.94, 134.35, 124.10, 146.22, 161.28, 151.72},
            {0.00, 21.34, 13.37, 21.03, 18.60, 10.68, 35.24,  27.24,  20.45,  27.09, 23.41,  28.14,  16.70,  28.91,  23.69,  20.73,  30.98,  18.53,  25.58,  22.22,  27.33,  22.78},
            {0.00, 27.88, 46.80, 34.87, 39.59, 51.30, 14.10,  8.63,   12.98,  11.19, 6.57,   7.19,   7.27,   7.75,   15.33,  17.10,  11.39,  6.06,   19.09,  19.56,  8.17,   12.95},
            {0.00, 43.41, 41.86, 43.46, 39.52, 37.67, 29.66,  20.64,  24.14,  24.60, 20.05,  16.14,  27.20,  18.07,  8.64,   10.28,  13.62,  22.47,  26.04,  30.79,  32.14,  16.72},
            {0.00, 24.75, 24.10, 20.72, 26.34, 32.54, 12.50,  10.07,  12.01,  10.10, 5.79,   9.59,   4.47,   9.07,   3.67,   2.72,   0.50,   4.17,   6.04,   4.68,   5.06,   6.86},
            {0.00, 24.50, 27.26, 32.70, 33.78, 68.57, 4.72,   8.28,   9.29,   5.69,  6.48,   9.81,   2.83,   4.72,   2.81,   4.30,   6.34,   4.39,   4.97,   2.31,   2.32,   4.48},
            {0.00, 65.12, 73.29, 63.89, 63.83, 66.10, 19.71,  47.48,  40.73,  46.87, 47.02,  44.40,  57.13,  48.13,  45.67,  51.01,  45.84,  42.78,  63.92,  42.27,  42.70,  42.62},
            {0.00, 62.22, 54.29, 45.62, 38.97, 46.68, 27.79,  40.41,  32.21,  38.42, 46.34,  31.20,  42.81,  35.87,  36.03,  43.75,  32.32,  29.64,  40.00,  38.39,  38.70,  34.80},
            {0.00, 40.07, 38.65, 28.91, 31.05, 26.75, 89.06,  93.31,  66.73,  91.34, 92.54,  85.32,  95.41,  104.47, 90.71,  93.89,  99.55,  100.21, 91.49,  83.91,  98.38,  98.44},
            {0.00, 16.82, 16.68, 14.72, 31.43, 16.29, 24.06,  32.17,  29.37,  30.31, 13.28,  17.98,  11.73,  9.99,   11.54,  7.14,   5.73,   9.61,   10.51,  5.04,   7.55,   15.28},
            {0.00, 22.92, 29.17, 36.67, 46.35, 19.61, 37.61,  28.35,  32.78,  36.96, 20.63,  21.80,  26.70,  23.77,  15.82,  11.39,  18.64,  12.12,  17.00,  13.91,  10.75,  20.42}
    };

    /*C-terminus of TM-helix (inside to outside)
    * matrix of size 20x21   (20 -> one for each amino acid)
    * columns -> position from TMS-center
    *  rows  -> amino acid order as ACDEFGHIKLMNPQRSTVWY
    */
    const std::vector<std::vector<double> > io_ct = {
            {0.00, 0.00,   0.00,   0.00,  0.00,   0.00,  0.00,  0.00,   0.00,   0.00,   0.00,  0.00,  0.00,   0.00,  0.00,  0.00,   0.00,  0.00,  0.00,  0.00,  0.00,  0.00},
            {0.00, 63.28,  52.81,  50.54, 63.44,  92.33, 59.76, 72.64,  53.93,  74.66,  59.40, 50.17, 69.36,  61.54, 58.57, 63.37,  43.05, 34.43, 38.80, 41.58, 33.81, 39.12},
            {0.00, 14.34,  18.58,  17.75, 20.55,  21.21, 22.34, 17.20,  20.00,  14.82,  15.42, 15.58, 12.79,  14.33, 8.60,  13.33,  8.64,  6.81,  4.18,  5.86,  10.56, 7.88},
            {0.00, 7.05,   3.70,   10.06, 7.13,   4.86,  4.09,  6.72,   5.22,   0.01,   2.73,  4.34,  5.02,   1.31,  5.42,  5.57,   5.45,  41.00, 34.19, 17.76, 20.07, 30.43},
            {0.00, 1.50,   2.40,   3.32,  5.25,   2.00,  7.67,  1.97,   1.12,   6.00,   2.08,  9.31,  5.29,   8.75,  7.25,  4.48,   3.00,  41.90, 23.46, 30.60, 35.86, 41.09},
            {0.00, 54.07,  55.00,  64.11, 64.24,  52.55, 67.27, 42.31,  56.71,  50.95,  48.51, 45.14, 50.35,  52.81, 49.59, 52.33,  56.05, 16.78, 28.18, 38.29, 23.23, 23.57},
            {0.00, 39.43,  50.94,  34.52, 28.76,  39.56, 54.75, 43.80,  47.97,  48.16,  49.21, 44.47, 38.64,  51.45, 44.28, 60.23,  54.22, 51.36, 42.18, 43.67, 51.28, 31.53},
            {0.00, 3.20,   2.00,   5.21,  1.33,   3.00,  2.01,  3.83,   5.08,   1.00,   10.71, 6.25,  7.08,   6.13,  5.19,  6.84,   4.47,  16.17, 20.51, 10.37, 15.06, 12.70},
            {0.00, 68.06,  72.65,  81.86, 63.95,  70.13, 69.43, 71.53,  59.72,  62.86,  64.91, 79.91, 72.53,  67.37, 52.06, 52.13,  56.01, 14.84, 22.45, 37.54, 41.04, 22.49},
            {0.00, 3.09,   2.03,   4.25,  2.50,   2.00,  2.00,  1.00,   5.50,   1.75,   4.05,  6.07,  2.08,   2.01,  8.76,  2.17,   1.74,  47.41, 24.96, 27.65, 29.21, 29.02},
            {0.00, 110.14, 116.79, 95.63, 111.47, 91.96, 91.40, 104.46, 109.05, 102.74, 89.57, 97.72, 101.93, 77.45, 83.22, 104.73, 97.20, 21.79, 38.10, 61.48, 60.05, 58.28},
            {0.00, 15.79,  20.06,  20.16, 22.06,  24.05, 18.67, 16.55,  29.78,  16.96,  26.34, 19.09, 19.63,  24.21, 27.79, 22.90,  27.95, 5.18,  19.61, 14.88, 16.11, 13.09},
            {0.00, 5.56,   12.97,  9.13,  8.29,   9.86,  10.16, 3.70,   8.81,   7.37,   15.87, 15.09, 13.59,  9.22,  14.55, 12.88,  10.54, 40.39, 29.97, 16.96, 20.43, 23.76},
            {0.00, 10.28,  8.62,   11.92, 8.68,   10.09, 17.73, 11.44,  10.09,  26.61,  24.69, 26.26, 25.70,  24.76, 16.82, 11.81,  16.75, 26.28, 40.16, 32.82, 34.66, 31.01},
            {0.00, 10.46,  6.83,   7.99,  7.97,   12.62, 9.60,  8.44,   11.34,  4.78,   7.53,  10.50, 8.78,   5.46,  14.54, 7.64,   8.39,  31.64, 16.72, 19.85, 18.01, 27.99},
            {0.00, 7.74,   1.18,   1.70,  1.00,   0.00,  1.00,  2.92,   4.70,   4.50,   2.83,  4.55,  6.36,   4.25,  4.92,  0.96,   4.25,  56.92, 29.33, 16.94, 20.04, 27.88},
            {0.00, 38.72,  35.53,  42.25, 30.52,  48.03, 35.37, 38.75,  36.21,  40.10,  37.22, 32.44, 27.01,  39.60, 43.12, 30.27,  24.46, 48.05, 46.70, 39.00, 40.88, 39.34},
            {0.00, 34.20,  29.28,  34.04, 22.96,  24.97, 32.05, 27.99,  34.15,  34.11,  32.77, 34.70, 34.40,  27.75, 32.16, 33.56,  22.78, 35.17, 47.15, 46.61, 38.77, 28.79},
            {0.00, 62.26,  65.14,  65.04, 77.78,  61.49, 63.14, 71.55,  59.53,  68.85,  60.19, 64.11, 53.98,  69.72, 59.25, 60.23,  69.48, 22.69, 36.68, 37.64, 38.52, 48.97},
            {0.00, 22.28,  19.40,  19.20, 24.09,  14.53, 7.47,  24.82,  10.90,  9.77,   14.19, 9.80,  14.63,  13.47, 14.58, 18.32,  31.70, 14.20, 15.51, 27.13, 17.94, 23.53},
            {0.00, 14.64,  13.17,  10.40, 18.10,  4.85,  14.18, 18.46,  20.27,  14.06,  21.85, 14.58, 20.92,  28.50, 39.41, 26.36,  43.96, 17.05, 31.25, 23.44, 24.56, 27.96}
    };


    /*C-terminus of TM-helix (outside to inside)
    * matrix of size 20x21   (20 -> one for each amino acid)
    * columns -> position from TMS-center
    *  rows  -> amino acid order as ACDEFGHIKLMNPQRSTVWY
    */
    const std::vector<std::vector<double> > oi_ct = {
            {0.00, 0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,  0.00,  0.00},
            {0.00, 86.65,  84.60,  90.94,  80.15,  83.81,  80.35,  71.75,  76.32,  73.97,  48.45,  67.42,  82.26,  76.92,  59.08,  73.61,  79.67,  36.05,  39.71,  37.90, 45.60, 53.44},
            {0.00, 17.72,  20.76,  22.50,  11.88,  20.60,  21.92,  23.81,  24.66,  18.71,  16.31,  26.03,  20.22,  16.16,  23.73,  13.33,  25.18,  14.56,  27.08,  18.30, 11.41, 9.19},
            {0.00, 8.05,   5.00,   2.67,   6.50,   2.55,   3.00,   3.62,   3.00,   4.75,   3.17,   4.38,   1.06,   1.19,   9.00,   2.30,   3.01,   44.44,  21.93,  23.79, 20.75, 21.40},
            {0.00, 8.30,   4.07,   4.00,   9.50,   4.47,   6.86,   10.41,  6.05,   5.45,   7.95,   4.25,   5.41,   4.63,   6.13,   8.12,   6.67,   38.09,  33.23,  28.35, 32.92, 33.56},
            {0.00, 48.08,  62.60,  79.57,  66.89,  74.47,  51.59,  58.30,  69.98,  56.20,  48.08,  56.56,  59.86,  63.30,  66.34,  80.32,  70.25,  18.28,  33.02,  41.78, 40.42, 33.08},
            {0.00, 73.33,  78.65,  70.65,  67.70,  51.99,  61.93,  55.24,  50.36,  60.09,  69.53,  54.15,  53.15,  46.91,  45.47,  39.45,  31.51,  33.31,  39.37,  46.21, 48.42, 41.12},
            {0.00, 2.13,   3.25,   2.83,   7.00,   3.80,   5.00,   1.90,   5.55,   0.90,   0.83,   3.02,   3.89,   6.31,   10.18,  3.73,   6.40,   30.10,  18.08,  16.91, 17.59, 23.70},
            {0.00, 86.84,  86.32,  68.55,  93.90,  96.53,  98.68,  86.54,  89.80,  97.95,  84.17,  92.48,  82.70,  90.42,  82.82,  89.23,  79.06,  19.63,  52.85,  40.22, 45.54, 42.09},
            {0.00, 7.58,   1.00,   4.83,   4.00,   1.52,   1.54,   3.20,   3.95,   4.56,   9.61,   8.24,   9.54,   3.33,   7.04,   7.83,   11.22,  123.75, 82.34,  73.05, 65.69, 61.53},
            {0.00, 113.83, 132.95, 136.12, 122.98, 145.77, 133.72, 147.32, 127.67, 135.38, 170.41, 154.98, 140.20, 123.78, 124.65, 132.06, 137.81, 36.24,  44.08,  57.08, 64.65, 65.04},
            {0.00, 19.42,  28.15,  18.64,  30.48,  21.19,  27.89,  25.58,  27.39,  23.54,  20.29,  31.02,  26.17,  20.63,  32.01,  19.36,  22.32,  8.53,   18.98,  24.77, 22.78, 14.33},
            {0.00, 8.00,   8.28,   18.50,  15.65,  17.58,  13.23,  12.00,  16.25,  30.23,  25.91,  10.63,  13.26,  8.81,   12.39,  10.35,  8.60,   42.92,  28.22,  33.40, 29.48, 24.51},
            {0.00, 20.17,  14.17,  13.07,  11.55,  13.43,  19.91,  35.84,  8.58,   25.43,  32.84,  20.47,  14.73,  16.21,  11.95,  11.77,  18.02,  19.08,  26.99,  24.31, 28.58, 25.44},
            {0.00, 5.92,   4.40,   3.87,   5.17,   4.00,   2.83,   5.33,   3.50,   9.83,   1.25,   9.85,   5.30,   5.58,   8.08,   10.63,  4.69,   24.58,  23.02,  22.97, 29.09, 26.71},
            {0.00, 5.82,   6.01,   2.00,   1.55,   4.85,   2.37,   5.52,   8.06,   1.25,   2.92,   14.04,  7.34,   11.14,  9.32,   14.74,  6.39,   114.83, 101.74, 80.13, 99.02, 87.31},
            {0.00, 50.00,  41.29,  48.43,  45.95,  39.89,  47.56,  55.70,  52.36,  45.27,  50.15,  32.37,  37.38,  32.47,  40.36,  45.00,  33.78,  46.26,  43.87,  43.49, 49.39, 52.52},
            {0.00, 34.78,  36.69,  28.19,  48.07,  40.58,  37.82,  27.98,  35.20,  33.20,  38.98,  33.20,  41.80,  28.74,  32.65,  29.40,  24.86,  31.61,  39.53,  39.47, 27.76, 34.82},
            {0.00, 93.40,  89.95,  96.58,  85.04,  88.04,  104.04, 96.75,  105.98, 83.43,  91.63,  85.66,  88.12,  109.83, 88.71,  81.05,  86.26,  25.09,  28.34,  34.63, 22.69, 43.13},
            {0.00, 12.71,  15.62,  12.82,  13.68,  8.69,   9.25,   7.09,   11.12,  8.57,   7.04,   10.55,  14.36,  19.21,  17.98,  23.35,  32.65,  11.38,  12.22,  12.92, 15.14, 17.89},
            {0.00, 34.41,  16.39,  16.41,  13.49,  17.42,  12.17,  7.75,   15.90,  22.93,  12.13,  22.37,  34.89,  56.09,  53.75,  46.03,  53.31,  22.91,  27.09,  41.97, 24.73, 27.84}
    };

    /* Normalized amino acid relative frequencies array
     */
    const std::vector<double> aa_freq_old = {
             0.00,
       0.083, 0.017, 0.053, 0.062, 0.039, 0.072, 0.022, 0.052, 0.057, 0.090,
       0.024, 0.044, 0.051, 0.040, 0.057, 0.069, 0.058, 0.066, 0.013, 0.032
    };

    /* Normalized amino acid relative frequencies array in TM-proteins
     */
    const std::vector<double> aa_freq_tm = {
             0.00,
            0.07441428, 0.02103942, 0.04459088, 0.05329261, 0.04564494,
            0.07045946, 0.02148279, 0.06092001, 0.0468635,  0.10051948,
            0.02460045, 0.04321853, 0.05013664, 0.0362897, 0.04713391,
            0.07238561, 0.06373683, 0.07077157, 0.01827555, 0.03422384
    };

    /* Normalized amino acid relative frequencies array (swissprot 25)
     */
    const std::vector<double> aa_freq = {
            0.00,
        0.0768, 0.0180, 0.0526, 0.0626, 0.0398, 0.0708, 0.0226, 0.0551, 0.0581, 0.0917,
        0.0235, 0.0443, 0.0505, 0.0403, 0.0526, 0.0707, 0.0584, 0.0652, 0.0131, 0.0322
    };

    /* Hydrophobicity (Eisenberg, Ann. Rev. Biochem. 53, 595 (1984))
     * Not used in tmpred
     */
    const std::vector<double> hphei = {
            0.00,
         0.25, -1.76, -0.64, -0.72, 0.04, -0.69, -0.62, 0.16, -0.4, 0.73,
         0.53, -1.1, 0.26, 0.61, -0.07, -0.26, -0.18, 0.37, 0.02, 0.54, 0
    };

    /* Hydrophobicity (Kyte and Doolittle, J. Mol. Biol. 157, 105 (1982))
     * Not used in tmpred
     */
    const std::vector<double> hypat = {
             0.00,
            1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5,
            3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2, 0
    };

    void scalematrix(std::vector<std::vector<double> > &d, std::vector<double> scale) {
        uint8_t ncol = n_col; // 21
        uint8_t nrow = n_row; // 20 one per each aminoacid

        for (uint8_t i = 1; i <= ncol; i++) {
            d[0][i] = 0;
            for (uint8_t j = 1; j <= nrow; j++)
                d[0][i] += d[j][i] / nrow;
        }

        double sum;
        for (uint8_t i = 1; i <= ncol; i++) {
            sum = 0;
            for (uint8_t j = 1; j <= nrow; j++)
                sum += d[j][i];
            for (uint8_t j = 1; j <= nrow; j++) {
                if (d[j][i] == 0)
                    d[j][i] = 1.0 / sum;  //can be improved?
                d[j][i] = d[j][i] / sum;
                //get the scale of the character (no need to be computed->should be j)
                d[j][i] = std::log(d[j][i] / scale[j]);
            }
            d[0][i] = 0;
            for (uint8_t j = 1; j <= nrow; j++)
                d[0][i] += d[j][i] / nrow;  //check if it used or not --> used when character is not part of rowcode
        }
    }

    void scaletables(std::vector<double> &scale,
                     std::vector<std::vector<double> > &m_io_ce,
                     std::vector<std::vector<double> > &m_oi_ce,
                     std::vector<std::vector<double> > &m_io_nt,
                     std::vector<std::vector<double> > &m_oi_nt,
                     std::vector<std::vector<double> > &m_io_ct,
                     std::vector<std::vector<double> > &m_oi_ct) {
        //first compute position 0 of the scale (aa_freq by default)
        scale[0] = 0;
        for (uint8_t i = 1; i < scale.size(); i++)
            scale[0] += scale[i] / n_row;
        m_io_ce = io_center;
        scalematrix(m_io_ce, scale);
        m_oi_ce = oi_center;
        scalematrix(m_oi_ce, scale);
        m_io_nt = io_nt;
        scalematrix(m_io_nt, scale);
        m_oi_nt = oi_nt;
        scalematrix(m_oi_nt, scale);
        m_io_ct = io_ct;
        scalematrix(m_io_ct, scale);
        m_oi_ct = oi_ct;
        scalematrix(m_oi_ct, scale);
    }



}


#endif //FDGPI_TABLES_H
