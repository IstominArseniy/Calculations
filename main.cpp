#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

// global constants
// ---------------------
double CHI_SCALE = 0.005;
double P_SCALE = 0.005;
double CHI_RANGE = 3.14159 / 2;
double P_RANGE = 2.0;
double P_DOT_SCALE = 3.0e-17;
int STEPS_NUMBER = 50;
int BIRTH_RATE = 1;
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

class Table{
private:
    int CHI_N;
    int P_N;
    double B12;
    double T_SCALE;
    table_t data_table;
    vector< vector< pair<double, double> > > dot_table;
    double find_T_scale(){
        double sum = 0.0;
        ofstream out("output.txt");
        for(int i = 1; i < P_N; i++){
            for(int j = 1; j < CHI_N; j++){
                out << dot_table[i][j].first << " ";
                if(death_line_check(i, j))
                    sum += max (2.0 * P_SCALE / dot_table[i][j].first, 2.0 * CHI_SCALE / dot_table[i][j].second);
            }
            out << '\n';
        }
        out.close();
        //cout << sum / P_N / CHI_N;
        return sum / P_N / CHI_N;
    }
    double birth_function(double P, double chi){
        if(MODEL == "BGI"){
            if(P < 0.5)
                return P * BGI_BIRTH_COEFFICIENT;
            else if(P < 1.0)
                return 0.5 * BGI_BIRTH_COEFFICIENT;
            else if( P < 1.5)
                return  1.5 * BGI_BIRTH_COEFFICIENT - P * BGI_BIRTH_COEFFICIENT;
            else
                return 0;
        }
        else if(MODEL == "MGD"){
            return sin(chi) * MGD_BIRTH_COEFFICIENT;
        }
    }

    bool death_line_check(int i, int j){
        /*
         return true if pulsar is alive
         */
        double P = (double)i * P_SCALE;
        double chi = (double)j * CHI_SCALE;
        //return true;
        if(pow(cos(chi), 0.4667) >= P * pow(A, 0.9333) * pow(B12, -0.5333))
            return true;
        return false;
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
        double d_P, d_chi;
        if(MODEL == "BGI"){
            d_P = dt * dot_table[i][j].first;
            d_chi = dt * dot_table[i][j].second;
        }
        else if(MODEL == "MGD"){
            ; //TODO do MGD
        }
        int new_i, new_j;
        new_i = (cr_P + d_P) / P_SCALE;
        new_j = (cr_chi + d_chi) / CHI_SCALE;
        //cout << new_i - i << " " << steps << " " << new_j - j << '\n';
        return make_pair(new_i, new_j);
    }

    void calculate_dot_table(){
        for(int i = 1; i < P_N; i++){
            for(int j = 1; j < CHI_N; j++){
                double cr_chi = (double)j * CHI_SCALE;
                double cr_P = (double)i * P_SCALE;
                if(MODEL == "BGI") {
                    dot_table[i][j].first = 10e-15 * B12 * B12 / cr_P *
                                            (Q(cr_P, cr_chi) * cos(cr_chi) * cos(cr_chi) + eps * pow(cr_P, -0.5));
                    dot_table[i][j].second = 10e-15 * Q(cr_P, cr_chi) * B12 * B12 / cr_P / cr_P * sin(cr_chi) * cos(cr_chi);
                }
                else{
                    return; // TODO MGD
                }
            }
        }
        return;
    }
    void update_table(int step){
        // Birth
        for(int i = 1; i < P_N; i++){
            for(int j = 1; j < CHI_N; j++) {
                data_table[i][j][(step + 1) % 2].number = 0.0;
                if(step % BIRTH_RATE == 0 and death_line_check(i, j)){
                    data_table[i][j][step % 2].number += birth_function((double)i * P_SCALE, (double)j * CHI_SCALE);
                }
            }
        }
        // Evolution
        for(int i = 1; i < P_N; i++){
            for(int j = 1; j < CHI_N; j++){
                pair<int, int> next_cell = find_next_cell(i, j, data_table[i][j][step % 2].steps_skipped + 1);
                //std::cout << i << " " << j << " " << next_cell.first << " " << next_cell.second <<  '\n';
                if(next_cell.first < P_N and next_cell.second < CHI_N and death_line_check(next_cell.first, next_cell.second)) {
                    data_table[next_cell.first][next_cell.second][(step + 1) % 2].number += data_table[i][j][step % 2].number;
                    if (next_cell.first != i or next_cell.second != j) {
                        data_table[next_cell.first][next_cell.second][(step + 1) % 2].steps_skipped = 0;
                    }
                    else{
                        data_table[next_cell.first][next_cell.second][(step + 1) % 2].steps_skipped = data_table[i][j][step % 2].steps_skipped + 1;
                    }
                }
                else{

                }
            }
        }
        return;
    }


public:
    Table(double B){
        CHI_N = CHI_RANGE / CHI_SCALE;
        P_N = P_RANGE / P_SCALE;
        B12 = B;
        data_table = table_t (P_N, vector< vector<cell> > (CHI_N, vector<cell>(2)));
        dot_table = vector< vector< pair<double, double> > > (P_N, vector <pair<double, double> >(CHI_N, make_pair(0.0, 0.0)));
        for(int i = 0; i < P_N; i++){
            for(int j = 0; j < CHI_N; j++){
                data_table[i][j][0].number = 0.0;
                data_table[i][j][1].number = 0.0;
                data_table[i][j][0].steps_skipped = 0.0;
                data_table[i][j][1].steps_skipped = 0.0;
            }
        }
        calculate_dot_table();
        T_SCALE = find_T_scale();
    }

    vector<vector<double> > get_P_chi_table(){
        vector<vector<double> > v(P_N - 1, vector<double> (CHI_N - 1, 0.0));
        for(int i = 1; i < P_N; i++){
            for(int j = 1; j < CHI_N; j++){
                v[i - 1][j - 1] = data_table[i][j][STEPS_NUMBER % 2].number;
            }
        }
        return v;
    }

    void calculate_P_chi_table(){
        for(int k = 0; k < STEPS_NUMBER; k++)
            update_table(k);
    }

    vector<vector<double> > get_P_P_dot_table(){
        double P_dot_RANGE = 1e-40;
        for(int i = 1; i < P_N; i++){
            for(int j = 1; j < CHI_N; j++){
                if(data_table[i][j][STEPS_NUMBER % 2].number != 0)
                    P_dot_RANGE = max(P_dot_RANGE, dot_table[i][j].first);
            }
        }
        int P_dot_N = P_dot_RANGE / P_DOT_SCALE;
        cout << P_dot_RANGE << " *\n";
        vector<vector<double> > v(P_N - 1, vector<double> (P_dot_N, 0.0));
        for(int i = 1; i < P_N; i++){
            for(int j = 1; j < CHI_N; j++){
                if(data_table[i][j][STEPS_NUMBER % 2].number != 0) {
                    int P_dot_ind = min(int(dot_table[i][j].first / P_dot_RANGE * P_dot_N), P_dot_N - 1);
                   // if (P_dot_ind > 100)
                      //  cout << P_dot_ind << " " << data_table[i][j][STEPS_NUMBER % 2].number <<  '\n';
                    v[i - 1][P_dot_ind] += data_table[i][j][STEPS_NUMBER % 2].number;
                }
            }
        }
        return v;
    }
};


void show_table(int N, int M, vector<vector<double > > const &v){
    ofstream out("data.txt");
    out << N << " " << M << '\n';
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            out << v[i][j] << " ";
        }
        out << '\n';
    }
    out.close();
    system("python plotter.py");
}


int main() {
    // magnetic fields
    double B12_1 = 1.0;
    double B12_2 = 3.0;
    double B12_3 = 5.0;
    double B12_4 = 7.0;
    double B12_5 = 10.0;
    //----------------
    Table table1(B12_1);
    Table table3(B12_3);
    //Table table5(B12_5);
    table1.calculate_P_chi_table();
    table3.calculate_P_chi_table();
    //table5.calculate_P_chi_table();
    vector<vector<double > > P_P_dot_table;
    vector<vector<double > > P_P_dot_table1 = table1.get_P_P_dot_table();
    vector<vector<double > > P_P_dot_table3 = table3.get_P_P_dot_table();
    //vector<vector<double > > P_CHI_table5 = table5.get_P_chi_table();
    int N = P_P_dot_table1.size();
    int M1 = P_P_dot_table1[0].size();
    int M3 = P_P_dot_table3[0].size();
    cout << M1 << " " << M3;
    int M = max(M1, M3);
    P_P_dot_table =  vector<vector<double > >(N, vector<double> (M, 0.0));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M1; j++){
            P_P_dot_table[i][j] += P_P_dot_table1[i][j];
        }
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M3; j++){
            P_P_dot_table[i][j] += P_P_dot_table3[i][j];
        }
    }

    show_table(N, M, P_P_dot_table);

    return 0;
}
