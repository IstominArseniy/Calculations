#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

// global constants
// ---------------------
double CHI_SCALE = 0.01;
double P_SCALE = 0.01;
double CHI_RANGE = 3.14159 / 2;
double P_RANGE = 1.0;
int STEPS_NUMBER = 10;
int BIRTH_RATE = 1;
double T_SCALE = 5.0e12; //TODO auto T scaling
double B12 = 1.0;
string MODEL = "BGI";
// BGI constants
//##############
double eps = 0.02;
double A = 1.0;
double D = 0.75;
int  BGI_BIRTH_COEFFICIENT = 10;
//##############
// MGD constants
//##############
int  MGD_BIRTH_COEFFICIENT = 1;
//##############
//-----------------------


struct cell{
    double number;
    int steps_skipped;
};


typedef vector< vector< vector<cell> > > table_t;



void reset_table(table_t& table, int P_N, int CHI_N){
    for(int i = 0; i < P_N; i++){
        for(int j = 0; j < CHI_N; j++){
            table[i][j][0].number = 0.0;
            table[i][j][1].number = 0.0;
            table[i][j][0].steps_skipped = 0.0;
            table[i][j][1].steps_skipped = 0.0;
        }
    }
    return;
}

double Q(double P, double chi){
    if(MODEL == "BGI"){
        double q = A * pow(P, 1.0714) * pow(B12, -0.5714) * pow(cos(chi), 2.0 * D - 2.0);
        if(q < 1.0)
            return A * pow(P, 1.0714) * pow(B12, -0.5714) * pow(cos(chi), 2.0 * D - 2.0); // question: dimension of Q?
        else
            return 1.0;
    }
    else if(MODEL == "MGD"){
        return 0.0; //TODO do MGD
    }
}

pair<int, int> find_next_cell(int i, int j, int steps){
    double cr_chi = (double)j * CHI_SCALE;
    double cr_P = (double)i * P_SCALE;
    double dt = (double)steps * T_SCALE;
    //std::cout << dt << '\n';
    double d_P, d_chi;
    if(MODEL == "BGI"){
        d_P = dt * 10e-15 * B12 * B12 / cr_P * (Q(cr_P, cr_chi) * cos(cr_chi) * cos(cr_chi) + eps * pow(cr_P, -0.5));
        d_chi = dt * 10e-15 * Q(cr_P, cr_chi) * B12 * B12 / cr_P / cr_P * sin(cr_chi) * cos(cr_chi);
    }
    else if(MODEL == "MGD"){
        ; //TODO do MGD
    }
    int new_i, new_j;
    //std::cout << d_P / P_SCALE << " " << d_chi / CHI_SCALE << '\n';
    new_i = (cr_P + d_P) / P_SCALE;
    new_j = (cr_chi + d_chi) / CHI_SCALE;
    return make_pair(new_i, new_j);
}

bool death_line_check(int i, int j){
    /*
     return true if pulsar is alive
     */
    double P = (double)i * P_SCALE;
    double chi = (double)j * CHI_SCALE;
    if(pow(cos(chi), 0.4667) >= P * pow(A, 0.9333) * pow(B12, -0.5333))
        return true;
    return false;
}
double birth_function(double P, double chi){
    if(MODEL == "BGI"){
        return P * BGI_BIRTH_COEFFICIENT;
    }
    else if(MODEL == "MGD"){
        return sin(chi) * MGD_BIRTH_COEFFICIENT;
    }
}
void update_table(table_t& table, int P_N, int CHI_N, int step){
    for(int i = 0; i < P_N; i++){
        for(int j = 0; j < CHI_N; j++) {
            table[i][j][(step + 1) % 2].number = 0;
            if(step % BIRTH_RATE == 0 and death_line_check(i, j)){
                table[i][j][step % 2].number += birth_function((double)i * P_SCALE, (double)j * CHI_SCALE);
            }
        }
    }
    for(int i = 0; i < P_N; i++){
        for(int j = 0; j < CHI_N; j++){
            pair<int, int> next_cell = find_next_cell(i, j, table[i][j][step % 2].steps_skipped + 1);
            // std::cout << i << " " << j << '\n';
            if(next_cell.first < P_N and next_cell.second < CHI_N and death_line_check(next_cell.first, next_cell.second)) {
                table[next_cell.first][next_cell.second][(step + 1) % 2].number += table[i][j][step % 2].number;
                if (next_cell.first != i or next_cell.second != j) {
                    table[next_cell.first][next_cell.second][(step + 1) % 2].steps_skipped = 0;
                }
                else{
                    table[next_cell.first][next_cell.second][(step + 1) % 2].steps_skipped =  table[i][j][(step + 1) % 2].steps_skipped + 1;
                }
            }
        }
    }
    return;
}




int main() {
    int CHI_N = CHI_RANGE / CHI_SCALE;
    int P_N = P_RANGE / P_SCALE;
    table_t table (P_N, vector< vector<cell> > (CHI_N, vector<cell>(2)));
    reset_table(table, P_N, CHI_N);
    // magnetic fields
    double B12_1 = 1.0;
    //----------------
    for(int i = 0; i < STEPS_NUMBER; i++){
        update_table(table, P_N, CHI_N, i);
    }

    ofstream out("data.txt");
    int Y_N, X_N;
    Y_N = P_N;
    X_N = CHI_N;
    out << P_N << " " << CHI_N << '\n';
    for(int i = 0; i < Y_N; i++){
        for(int j = 0; j < X_N; j++){
            out << table[i][j][STEPS_NUMBER % 2].number << " ";
        }
        out << "\n";
    }
    out.close();
    system("python plotter.py");
    return 0;
}
//
// Created by Arseniy on 30.11.2020.
//

