//**********************************************************************
// Gauss Konrod Quadrature weights and abscisas
//********************************************************************** 
void gauss_konrod_w(double GKQ1[2*nw1 + 1][3], 
                    double GKQ2[2*nw2 + 1][3],
                    double GKQ3[2*nw3 + 1][3],
                    double GKQ4[2*nw4 + 1][3]){

    int const N1 = nw1;
    int const N2 = nw2;
    int const N3 = nw3;
    int const N4 = nw4;

    auto XGK1 = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::abscissa();    
    auto WGK1 = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::weights();      
    auto WG1 =  boost::math::quadrature::gauss<double, N1>::weights();

    auto XGK2 = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::abscissa();
    auto WGK2 = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::weights();
    auto WG2 =  boost::math::quadrature::gauss<double, N2>::weights();

    auto XGK3 = boost::math::quadrature::gauss_kronrod<double, 2*N3+1>::abscissa();    
    auto WGK3 = boost::math::quadrature::gauss_kronrod<double, 2*N3+1>::weights();      
    auto WG3 =  boost::math::quadrature::gauss<double, N3>::weights();

    auto XGK4 = boost::math::quadrature::gauss_kronrod<double, 2*N4+1>::abscissa();
    auto WGK4 = boost::math::quadrature::gauss_kronrod<double, 2*N4+1>::weights();
    auto WG4 =  boost::math::quadrature::gauss<double, N4>::weights();

    double WGL1[2*N1+1], WGL2[2*N2+1], WGL3[2*N3+1], WGL4[2*N4+1];


// Changes the n order array of Gaussian quadrature to a 2n+1 array

    for(int i = 0; i < XGK1.size(); ++i){
        if (i % 2 == 0){
            WGL1[i] = WG1[i/2];}
        else WGL1[i] = 0.;
    }

    for(int i = 0; i < XGK2.size(); ++i){
        if (i % 2 == 0){
            WGL2[i] = WG2[i/2];}
        else WGL2[i] = 0.;
    }

    for(int i = 0; i < XGK3.size(); ++i){
        if (i % 2 == 0){
            WGL3[i] = WG3[i/2];}
        else WGL3[i] = 0.;
    }

    for(int i = 0; i < XGK4.size(); ++i){
        if (i % 2 == 0){
            WGL4[i] = WG4[i/2];}
        else WGL4[i] = 0.;
    }

// writtes [0] abscisas, [1] Konrod weights, [2] Gauss weigths

    for(int i = 0; i < 2*XGK1.size() - 1; ++i){
        if (i <= N1){ 
            GKQ1[i][0] = -1.*XGK1[N1 - i];
            GKQ1[i][1] = WGK1[N1 - i];                
            GKQ1[i][2] = WGL1[N1 - i];
        }
        else{
            GKQ1[i][0] = XGK1[-N1 + i];
            GKQ1[i][1] = WGK1[-N1 + i];
            GKQ1[i][2] = WGL1[-N1 + i];
        }
    }

    for(int i = 0; i < 2*XGK2.size() - 1; ++i){
        if (i <= N2){ 
            GKQ2[i][0] = -1.*XGK2[N2 - i];
            GKQ2[i][1] = WGK2[N2 - i];
            GKQ2[i][2] = WGL2[N2 - i];
        }
        else{
            GKQ2[i][0] = XGK2[-N2 + i];
            GKQ2[i][1] = WGK2[-N2 + i];
            GKQ2[i][2] = WGL2[-N2 + i];
        }
    }

    for(int i = 0; i < 2*XGK3.size() - 1; ++i){
        if (i <= N3){ 
            GKQ3[i][0] = -1.*XGK3[N3 - i];
            GKQ3[i][1] = WGK3[N3 - i];                
            GKQ3[i][2] = WGL3[N3 - i];
        }
        else{
            GKQ3[i][0] = XGK3[-N3 + i];
            GKQ3[i][1] = WGK3[-N3 + i];
            GKQ3[i][2] = WGL3[-N3 + i];
        }
    }

    for(int i = 0; i < 2*XGK4.size() - 1; ++i){
        if (i <= N4){ 
            GKQ4[i][0] = -1.*XGK4[N4 - i];
            GKQ4[i][1] = WGK4[N4 - i];
            GKQ4[i][2] = WGL4[N4 - i];
        }
        else{
            GKQ4[i][0] = XGK4[-N4 + i];
            GKQ4[i][1] = WGK4[-N4 + i];
            GKQ4[i][2] = WGL4[-N4 + i];
        }
    }

}

//******************************************************************************************************
//******************************************************************************************************

void Omegas(double xi[2*NN + 4], double xk[2*NN + 4], double xg[2*NN + 4]){

double GKQw1[2*nw1 + 1][3], GKQw2[2*nw2 + 1][3];
double GKQw3[2*nw3 + 1][3], GKQw4[2*nw4 + 1][3];
gauss_konrod_w(GKQw1,GKQw2,GKQw3,GKQw4);

double w;

int NN1 = 2*nw1 + 1;
int NN2 = 2*nw2 + 1;
int NN3 = 2*nw3 + 1;
int NN4 = 2*nw4 + 1;

for(int i=0; i < NN1;  ++i){
xi[i] = ((w2-w1)/2.)*GKQw1[i][0] + ((w2+w1)/2.);;
xk[i] = ((w2-w1)/2.)*GKQw1[i][1];
xg[i] = ((w2-w1)/2.)*GKQw1[i][2];
}

for(int i=0; i < NN2;  ++i){
xi[i + NN1] = ((w3-w2)/2.)*GKQw2[i][0] + ((w3+w2)/2.);;
xk[i + NN1] = ((w3-w2)/2.)*GKQw2[i][1];
xg[i + NN1] = ((w3-w2)/2.)*GKQw2[i][2];
}

for(int i=0; i < NN3;  ++i){
xi[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][0] + ((w4+w3)/2.);;
xk[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][1];
xg[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][2];
}

for(int i=0; i < NN4;  ++i){
xi[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][0] + ((w5+w4)/2.);;
xk[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][1];
xg[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][2];
}

}