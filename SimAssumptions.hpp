#pragma once
const static double MIN_AGE = 0.001;
const static double PLANET_FORM_AGE = 0.006;
const static double PRECISION = 1e-5;
const static double MAX_END_AGE = 10;
const static double SPIN_THRES = 4*M_PI; //0.5 days
const static double MAIN_SEQ_START = 0.2;
const static double WIND_SAT_FREQ = 2.454; //10x solar rotation rate
const static double Q_TRANS_WIDTH = 1e-3;

#define Q1 1e6
#define Q2 5e6
#define Q3 1e7
#define N_Q 3
#define N_ROT 3
#define N_SMASS 10
#define N_PMASS 50
#define N_P0 100

const double MIN_SMASS = 0.5, MAX_SMASS = 1.1;
const double MIN_PMASS = log(0.05), MAX_PMASS = log(25);
const double MIN_P0 = 0.1, MAX_P0 = 5.9;
