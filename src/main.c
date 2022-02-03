#include "math.h"
#include "unistd.h"
#include "stdio.h"
#include "basic_linear.h"

// Debug print
#define DBG 0
#define DEBUGPRT(x,a,b) if(DBG){printf(x,a,b);}

// Some Constant
#define PI 3.14159
#define MAX 65536
// Some config
#define CONV 1e-8
#define MAXMAINITER 10000
#define MAXMEQNITER 200
#define MAXPEQNITER 200
#define MAXPCORITER 200

// Define BC type
#define WALL   0 // fixed u and v
#define INLET  1 // floating u, fixed v
#define OUTLET 2 // floating u and v

// Define uorv value
#define ISU 1
#define ISV 0

// Define hat value
#define HAT 1
#define NOTHAT 0 // not hat

typedef struct BoundaryCondition{
    int type;
    double u0;
    double v0;
    double P0;
}BC;

typedef struct RectFlowField{
    double dx, dy, dt;
    double rho, miu;
    int M, N, L, H;
    double** U,** V,** Uhat,** Vhat,** P,** uap,** vap; // ap is the coeff. in momentum eqn.
    BC Nbc, Ebc, Sbc, Wbc; 
}RFF;

int CNT;

void setVelocityConstant(RFF* field);
void initialize_flow_field(RFF* field, int M, int N, double L, double H);
double UorUhat(RFF* field,int i,int j,int hat);
double VorVhat(RFF* field,int i,int j,int hat);
double get_u(RFF* field,int i,int j, int hat);
double get_v(RFF* field,int i,int j, int hat);
double theA(double P);
void momentum_equation_coeff(RFF* field, int i, int j, double* a, int hat, int uorv);
void pseudoU(RFF* field);
void pseudoV(RFF* field);
double get_ap(RFF* field, int i, int j, int uorv);
double check_P_boundary(RFF* field, int i, int j);
void pressure_equation_coeff(RFF* field, int i, int j, double* a);
void SolvePressureEquation(RFF* field);
double friction(RFF* field, int i, int j, int uorv);
void SolveMomentumEquation(RFF* field);
void SolvePressureCorrection(RFF* field, double** Pcorr);
double VelocityCorrection(RFF* field, double** Pcorr);

void couette_flow(void);
void cavity_flow(void);

int main(int argc, char **argv){
    if(argv[1][0]=='a'){
        couette_flow();
    }
    else{
        cavity_flow();
    }
}

void couette_flow(){
    int count = 0;
    int period = 1000;
    double diff=MAX;
    // Setting Boundary Condition
    BC Nth, Sth, Est, Wst;
    Nth.type = WALL;
    Nth.v0 = 0;
    Nth.u0 = 0.01;
    Sth.type = WALL;
    Sth.v0 = 0;
    Sth.u0 = 0;
    Wst.type = INLET;
    Wst.P0 = 0;
    Est.type = OUTLET;
    Est.v0 = 0;
    Est.P0 = 0;
    // Initialize flow field
    RFF field;
    field.rho = 997; // Kg/m3 for water
    field.miu = 0.001; // Pa*s fir water
    field.dt = 0.1;

    field.Nbc = Nth;
    field.Sbc = Sth;
    field.Ebc = Est;
    field.Wbc = Wst;
    initialize_flow_field(&field,10,10,1,0.025);
    // Initialize some temperal stuff
    double** Pcorr = dmatrix(field.M,field.N);
    //printf("Initialize Done \r\n");
    //printf("North u: %lf\r\n",field.Nbc.u0);
    printf("%d %d\n",field.M,field.N);
    while((diff>CONV)&&(count<MAXMAINITER)){
        CNT = count;
        pseudoU(&field);
        pseudoV(&field);
        SolvePressureEquation(&field);
        /*if(count==15000){
            printdmat(field.P,field.N,field.M);
        }*/
        SolveMomentumEquation(&field);
        makezero(Pcorr,field.M,field.N);
        SolvePressureCorrection(&field, Pcorr);
        diff = VelocityCorrection(&field, Pcorr);
        count++;
        /*if((count%period)==0){
            printf("Main Iteration: %d, diff: %lf\r\n",count,diff);
        }*//*
        if(count==5000){
            printdmat(Pcorr,field.N,field.M);
        }*/
    }
    printdmat(field.U,field.N+1,field.M);
    printdmat(field.V,field.N+2,field.M+1);
    printdmat(field.P,field.N,field.M);
}

void cavity_flow(){
    int count = 0;
    int period = 100;
    double diff=MAX;
    // Setting Boundary Condition
    BC Nth, Sth, Est, Wst;
    Nth.type = WALL;
    Nth.v0 = 0;
    Nth.u0 = 0.1;
    Sth.type = WALL;
    Sth.v0 = 0;
    Sth.u0 = 0;
    Wst.type = WALL;
    Wst.u0 = 0;
    Wst.v0 = 0;
    Est.type = WALL;
    Est.u0 = 0;
    Est.v0 = 0;
    // Initialize flow field
    RFF field;
    field.rho = 1.23; // Kg/m3 for air
    field.miu = 1.778e-5; // Pa*s (kg/m/s) for air
    field.dt = 0.01;

    field.Nbc = Nth;
    field.Sbc = Sth;
    field.Ebc = Est;
    field.Wbc = Wst;
    double box = 0.01;
    initialize_flow_field(&field,10,10,box,box);
    // Initialize some temperal stuff
    double** Pcorr = dmatrix(field.M,field.N);
    //printf("Initialize Done \r\n");
    //printf("Re: %lf\r\n",field.rho*Nth.u0*box/field.miu);
    while((diff>CONV)&&(count<MAXMAINITER)){
        CNT = count;
        pseudoU(&field);
        pseudoV(&field);
        SolvePressureEquation(&field);
        SolveMomentumEquation(&field);
        makezero(Pcorr,field.M,field.N);
        SolvePressureCorrection(&field, Pcorr);
        diff = VelocityCorrection(&field, Pcorr);
        count++;
        /*if((count%period)==0){
            printf("Main Iteration: %d, diff: %lf\r\n",count,diff);
        }*/
    }
    printf("%d %d\n",field.M,field.N);
    printdmat(field.U,field.N+1,field.M);
    printdmat(field.V,field.N+2,field.M+1);
    printdmat(field.P,field.N,field.M);

}

void setVelocityConstant(RFF* field){
    // If Wall conditions are set, we need to set velocities back to constants at walls.
    // East
    if(field->Ebc.type==WALL){
        for(int i=0;i<field->M;i++){
            field->U[field->N][i] = field->Ebc.u0;
            field->Uhat[field->N][i] = field->Ebc.u0;
        }
    }
    if(field->Wbc.type==WALL){
        for(int i=0;i<field->M;i++){
            field->U[0][i] = field->Wbc.u0;
            field->Uhat[0][i] = field->Wbc.u0;
        }
    }
    if(field->Nbc.type==WALL){
        for(int i=0;i<field->N+2;i++){
            field->V[i][field->M] = field->Nbc.v0;
            field->Vhat[i][field->M] = field->Nbc.v0;
        }
    }
    if(field->Sbc.type==WALL){
        for(int i=0;i<field->N+2;i++){
            field->V[i][0] = field->Sbc.v0;
            field->Vhat[i][0] = field->Sbc.v0;
        }
    }
}

void initialize_flow_field(RFF* field, int M, int N, double L, double H){
    field->M = M;
    field->N = N;
    field->L = L;
    field->H = H;
    field->dx = L/N;
    field->dy = H/M;
    field->P = dzero(N, M);
    field->U = dzero(N+1, M);
    field->V = dzero(N+2,M+1);
    field->Uhat = dzero(N+1, M);
    field->Vhat = dzero(N+2, M+1);
    field->uap = dzero(N+1, M);
    field->vap = dzero(N+2,M+1);
    setVelocityConstant(field);
}

double UorUhat(RFF* field,int i,int j,int hat){
    if(hat){return field->Uhat[i][j];}
    else{return field->U[i][j];}
}

double VorVhat(RFF* field,int i,int j,int hat){
    if(hat){return field->Vhat[i][j];}
    else{return field->V[i][j];}
}

double get_u(RFF* field,int i,int j, int hat){
    if(i>field->N){ // Exceed East bound ((N+1)
        if(field->Ebc.type==WALL){
            return 0; // Don't consider ae and ae*ue
        }
        else if(field->Ebc.type==INLET){
            return get_u(field,i-1,j, hat); // Zero order extrapolate
        }
        else if(field->Ebc.type==OUTLET){
            return get_u(field,i-1,j, hat); // Zero order extrapolate
        }
        else{
            printf("Error East Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(i<0){ // Exceed West Bound
        if(field->Wbc.type==WALL){
            return 0; // Don't consider aw and aw*uw
        }
        else if(field->Wbc.type==INLET){
            return get_u(field,i+1,j, hat); // Zero order extrapolate
        }
        else if(field->Wbc.type==OUTLET){
            return get_u(field,i+1,j, hat); // Zero order extrapolate
        }
        else{
            printf("Error West Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j>=field->M){ // Exceed North Bound (M)
        if(field->Nbc.type==WALL){
            return 0; // Don't consider an and an*un
        }
        else if(field->Nbc.type==INLET){
            return field->Nbc.u0; // Zero order extrapolate
        }
        else if(field->Nbc.type==OUTLET){
            return get_u(field,i,j-1, hat); // Zero order extrapolate
        }
        else{
            printf("Error North Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j<0){
        if(field->Sbc.type==WALL){ // Exceed South Bound
            return 0; // Don't consider as and as*us
        }
        else if(field->Sbc.type==INLET){
            return field->Sbc.u0; // Zero order extrapolate
        }
        else if(field->Sbc.type==OUTLET){
            return get_u(field,i,j+1, hat); // Zero order extrapolate
        }
        else{
            printf("Error South Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }

    }
    else if(i==field->N){ // On East bound ((N+1)
        if(field->Ebc.type==WALL){
            return field->Ebc.u0; // Don't consider ae and ae*ue
        }
        else if(field->Ebc.type==INLET){
            return UorUhat(field,i,j,hat);
        }
        else if(field->Ebc.type==OUTLET){
            return UorUhat(field,i,j,hat);
        }
        else{
            printf("Error East Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(i==0){ // On West Bound
        if(field->Wbc.type==WALL){
            return field->Wbc.u0;
        }
        else if(field->Wbc.type==INLET){
            return UorUhat(field,i,j,hat);
        }
        else if(field->Wbc.type==OUTLET){
            return UorUhat(field,i,j,hat);
        }
        else{
            printf("Error West Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else{ // Interior Node
        return UorUhat(field,i,j,hat);
    }
}

double get_v(RFF* field,int i,int j, int hat){
    if(i>=field->N+1){ // Exceed East bound (N+1)
        if(field->Ebc.type==WALL){
            return 0; // Don't consider ae and ae*ue
        }
        else if(field->Ebc.type==INLET){
            return field->Ebc.v0; // Zero order extrapolate
        }
        else if(field->Ebc.type==OUTLET){
            return get_v(field,i-1,j, hat); // Zero order extrapolate
        }
        else{
            printf("Error East Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(i<1){ // Exceed West Bound
        if(field->Wbc.type==WALL){
            return 0; // Don't consider aw and aw*uw
        }
        else if(field->Wbc.type==INLET){
            return field->Wbc.v0; // Zero order extrapolate
        }
        else if(field->Wbc.type==OUTLET){
            return get_v(field,i+1,j, hat); // Zero order extrapolate
        }
        else{
            printf("Error West Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j>=field->M+1){ // Exceed North Bound (M+1)
        if(field->Nbc.type==WALL){
            return 0; // Don't consider an and an*un
        }
        else if(field->Nbc.type==INLET){
            return get_v(field,i,j-1, hat); // Zero order extrapolate
        }
        else if(field->Nbc.type==OUTLET){
            return get_v(field,i,j-1, hat); // Zero order extrapolate
        }
        else{
            printf("Error North Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j<0){ // Exceed South Bound
        if(field->Sbc.type==WALL){
            return 0; // Don't consider as and as*us
        }
        else if(field->Sbc.type==INLET){
            return get_v(field,i,j+1, hat); // Zero order extrapolate
        }
        else if(field->Sbc.type==OUTLET){
            return get_v(field,i,j+1, hat); // Zero order extrapolate
        }
        else{
            printf("Error South Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j==field->M){ // On North Bound (M+1)
        if(field->Nbc.type==WALL){
            return field->Nbc.v0; // Don't consider an and an*un
        }
        else if(field->Nbc.type==INLET){
            return VorVhat(field,i,j,hat);
        }
        else if(field->Nbc.type==OUTLET){
            return VorVhat(field,i,j,hat);
        }
        else{
            printf("Error North Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j==0){ // On South Bound
        if(field->Sbc.type==WALL){
            return 0; // Don't consider as and as*us
        }
        else if(field->Sbc.type==INLET){
            return VorVhat(field,i,j,hat);
        }
        else if(field->Sbc.type==OUTLET){
            return VorVhat(field,i,j,hat);
        }
        else{
            printf("Error South Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else{ // Interior Node
        return VorVhat(field,i,j,hat);
    }
}

double get_P(RFF* field, int i, int j, double** Pfield){
    if(i>=field->N){ // Exceed East bound (N)
        if(field->Ebc.type==WALL){
            return get_P(field,i-1,j,Pfield); // Don't consider ae and ae*ue
        }
        else if(field->Ebc.type==INLET){
            return 0; // Zero order extrapolate
        }
        else if(field->Ebc.type==OUTLET){
            return 0; // Zero order extrapolate
        }
        else{
            printf("Error East Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(i<0){ // Exceed West Bound
        if(field->Wbc.type==WALL){
            return get_P(field,i+1,j,Pfield); // Don't consider aw and aw*uw
        }
        else if(field->Wbc.type==INLET){
            return 0; // Zero order extrapolate
        }
        else if(field->Wbc.type==OUTLET){
            return 0; // Zero order extrapolate
        }
        else{
            printf("Error West Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j>=field->M){ // Exceed North Bound (M)
        if(field->Nbc.type==WALL){
            return get_P(field,i,j-1,Pfield); // Don't consider an and an*un
        }
        else if(field->Nbc.type==INLET){
            return 0; // Zero order extrapolate
        }
        else if(field->Nbc.type==OUTLET){
            return 0; // Zero order extrapolate
        }
        else{
            printf("Error North Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else if(j<0){ // Exceed South Bound
        if(field->Sbc.type==WALL){
            return get_P(field,i,j+1,Pfield); // Don't consider as and as*us
        }
        else if(field->Sbc.type==INLET){
            return 0; // Zero order extrapolate
        }
        else if(field->Sbc.type==OUTLET){
            return 0; // Zero order extrapolate
        }
        else{
            printf("Error South Boundary Type: (i,j)=(%d,%d)\r\n",i,j);
            return 0;
        }
    }
    else{ // Interior Node
        return Pfield[i][j];
    }
}

double theA(double P){
    return fmax(0,pow((1-0.1*abs(P)),5));
}

void momentum_equation_coeff(RFF* field, int i, int j, double* a,int hat, int uorv){
    double ap0 = field->rho*field->dx*field->dy/field->dt;
    if(uorv==ISU){
        double ue = 0.5*(get_u(field,i+1,j,hat)+get_u(field,i,j,hat));
        double uw = 0.5*(get_u(field,i,j,hat)+get_u(field,i-1,j,hat));
        double vn = 0.5*(get_v(field,i,j+1,hat)+get_v(field,i+1,j+1,hat));
        double vs = 0.5*(get_v(field,i,j,hat)+get_v(field,i+1,j,hat));

        double Fe = field->rho*ue*field->dy;
        double De = field->miu*field->dy/field->dx;
        double Pe = Fe/De;
        double ae = De*theA(Pe)+fmax(-Fe, 0);

        double Fw = field->rho*uw*field->dy;
        double Dw = field->miu*field->dy/field->dx;
        double Pw = Fw/Dw;
        double aw = Dw*theA(Pw)+fmax(Fw, 0);

        double Fn = field->rho*vn*field->dx;
        double Dn = field->miu*field->dx/field->dy;
        double Pn = Fn/Dn;
        double an = Dn*theA(Pn)+fmax(-Fn, 0);

        double Fs = field->rho*vs*field->dx;
        double Ds = field->miu*field->dx/field->dy;
        double Ps = Fs/Ds;
        double as = Ds*theA(Ps)+fmax(Fs, 0);
        // At wall condition, ax is 0 if it is parallel to wall
        if((j==field->M-1)&&(field->Nbc.type==WALL)){
            an = 0;
        }
        else if((j==0)&&(field->Sbc.type==WALL)){
            as = 0;
        }
        field->uap[i][j] = ae + aw + an + as + ap0;
        a[4] = field->uap[i][j];
        
        a[0] = ae;
        a[1] = aw;
        a[2] = an;
        a[3] = as;
        a[5] = ap0;
    }
    else if(uorv==ISV){
        double ue = 0.5*(get_u(field,i,j,hat)+get_u(field,i,j-1,hat));
        double uw = 0.5*(get_u(field,i-1,j,hat)+get_u(field,i-1,j-1,hat));
        double vn = 0.5*(get_v(field,i,j+1,hat)+get_v(field,i,j,hat));
        double vs = 0.5*(get_v(field,i,j,hat)+get_v(field,i,j-1,hat));

        double Fe = field->rho*ue*field->dy;
        double De = field->miu*field->dy/field->dx;
        double Pe = Fe/De;
        double ae = De*theA(Pe)+fmax(-Fe, 0);

        double Fw = field->rho*uw*field->dy;
        double Dw = field->miu*field->dy/field->dx;
        double Pw = Fw/Dw;
        double aw = Dw*theA(Pw)+fmax(Fw, 0);

        double Fn = field->rho*vn*field->dx;
        double Dn = field->miu*field->dx/field->dy;
        double Pn = Fn/Dn;
        double an = Dn*theA(Pn)+fmax(-Fn, 0);

        double Fs = field->rho*vs*field->dx;
        double Ds = field->miu*field->dx/field->dy;
        double Ps = Fs/Ds;
        double as = Ds*theA(Ps)+fmax(Fs, 0);
        // At wall condition, ax is 0 if it is parallel to wall
        if((i==field->M+1)&&(field->Ebc.type==WALL)){
            ae = 0;
        }
        else if((i==0)&&(field->Wbc.type==WALL)){
            aw = 0;
        }
        field->vap[i][j] = ae + aw + an + as + ap0;
        a[4] = field->vap[i][j];
        
        a[0] = ae;
        a[1] = aw;
        a[2] = an;
        a[3] = as;
        a[5] = ap0;
    }
}

void pseudoU(RFF* field){
    double a[6];
    for(int i=0;i<field->N+1;i++){
        for(int j=0;j<field->M;j++){
            momentum_equation_coeff(field, i, j, a, NOTHAT,ISU);
            double b = a[5]*get_u(field,i,j,NOTHAT); // b = phi0*a0
            field->Uhat[i][j] = (a[0]*get_u(field,i+1,j,NOTHAT)
                                +a[1]*get_u(field,i-1,j,NOTHAT)
                                +a[2]*get_u(field,i,j+1,NOTHAT)
                                +a[3]*get_u(field,i,j-1,NOTHAT)
                                +b)/a[4]
                                ;
        }
    }
    setVelocityConstant(field);
}

void pseudoV(RFF* field){
    double a[6];
    for(int i=0;i<field->N+2;i++){
        for(int j=0;j<field->M+1;j++){
            momentum_equation_coeff(field, i, j, a, NOTHAT,ISV);
            double b = a[5]*get_v(field,i,j,NOTHAT); // b = phi0*a0
            field->Vhat[i][j] = (a[0]*get_v(field,i+1,j,NOTHAT)
                                +a[1]*get_v(field,i-1,j,NOTHAT)
                                +a[2]*get_v(field,i,j+1,NOTHAT)
                                +a[3]*get_v(field,i,j-1,NOTHAT)
                                +b)/a[4]
                                ;
        }
    }
    setVelocityConstant(field);
}

double get_ap(RFF* field, int i, int j, int uorv){
    if((i<field->N)&&(i>=0)&&(j<field->M)&&(j>=0)){
        if(uorv == ISU){
            return field->uap[i][j];
        }
        else if(uorv == ISV){
            return field->vap[i][j];
        }
    }
    else{
        return 1; // It is at the denominator, so return 1
    }
}

double check_P_boundary(RFF* field, int i, int j){
    if((i<field->N)&&(i>=0)&&(j<field->M)&&(j>=0)){
        return 1;
    }
    else{
        return 0;
    }
}

void pressure_equation_coeff(RFF* field, int i, int j, double* a){
    double ae = check_P_boundary(field,i+1,j)*field->rho*field->dy*field->dy/get_ap(field,i+1,j,ISU);
    double aw = check_P_boundary(field,i-1,j)*field->rho*field->dy*field->dy/get_ap(field,i-1,j,ISU);
    double an = check_P_boundary(field,i,j+1)*field->rho*field->dx*field->dx/get_ap(field,i,j+1,ISV);
    double as = check_P_boundary(field,i,j-1)*field->rho*field->dx*field->dx/get_ap(field,i,j-1,ISV);
    double ap = ae+aw+an+as;
    a[0] = ae;
    a[1] = aw;
    a[2] = an;
    a[3] = as;
    a[4] = ap;
}

void SolvePressureEquation(RFF* field){
    int count = 0;
    double diff = MAX;
    makezero(field->P, field->N, field->M);
    while((count<=MAXPEQNITER)&&(diff>CONV)){
        diff = 0;
        for(int i=0;i<field->N;i++){
            for(int j=0;j<field->M;j++){
                double a[5];
                pressure_equation_coeff(field,i,j,a);
                double b = field->dy*(field->rho*get_u(field,i  ,j,HAT)
                                    - field->rho*get_u(field,i+1,j,HAT))
                         + field->dx*(field->rho*get_v(field,i+1,j,HAT)
                                    - field->rho*get_v(field,i+1,j+1,HAT));
                //printf("b in PEQN: %lf\r\n",b);
                double newPs = (a[0]*get_P(field,i+1,j,field->P) 
                             + a[1]*get_P(field,i-1,j,field->P)
                             + a[2]*get_P(field,i,j+1,field->P)
                             + a[3]*get_P(field,i,j-1,field->P)
                             + b)/a[4];
                diff += pow(newPs-field->P[i][j],2);
                field->P[i][j] = newPs;
            }
        }
        diff = sqrt(diff);
        DEBUGPRT("Solving Pressure Eqn. Iter: %d, Diff: %lf\r\n", count, diff)
        count++;
    }
}

double friction(RFF* field, int i, int j, int uorv){
    // Produce friction only at walls
    // when uorv == 1. it means u
    if((uorv==ISU)&&(j==field->M-1)){
        return field->dx*field->miu*2*(field->Nbc.u0-get_u(field,i,j,HAT))/field->dy; // dx*miu*(du/dy)
    }
    else if((uorv==ISU)&&(j==0)){
        return field->dx*field->miu*2*(field->Sbc.u0-get_u(field,i,j,HAT))/field->dy; // dx*miu*(du/dy)
    }
    else if((uorv==ISV)&&(i==field->N+1)){
        return field->dy*field->miu*2*(field->Ebc.v0-get_v(field,i,j,HAT))/field->dx; // dy*miu*(dv/dx)
    }
    else if((uorv==ISV)&&(i==0)){
        return field->dy*field->miu*2*(field->Wbc.v0-get_v(field,i,j,HAT))/field->dx; // dy*miu*(dv/dx)
    }
}

void SolveMomentumEquation(RFF* field){
    int count = 0;
    double diff = MAX;
    makezero(field->Uhat,field->M,field->N+1);
    makezero(field->Uhat,field->M+1,field->N+2);
    setVelocityConstant(field);
    // U first
    while((count<=MAXMEQNITER)&&(diff>CONV)){
        diff = 0;
        for(int i=0;i<field->N+1;i++){ 
            for(int j=0;j<field->M;j++){
                double a[6];
                momentum_equation_coeff(field,i,j,a,NOTHAT,ISU);
                double b = a[5]*get_u(field,i,j,NOTHAT);
                double newUs = (a[0]*get_u(field,i+1,j,HAT) 
                             + a[1]*get_u(field,i-1,j,HAT)
                             + a[2]*get_u(field,i,j+1,HAT)
                             + a[3]*get_u(field,i,j-1,HAT)
                             + b + field->dy*(get_P(field,i-1,j,field->P)-get_P(field,i,j,field->P))
                             + friction(field,i,j,ISU))/a[4];
                diff += pow(newUs-field->Uhat[i][j],2);
                field->Uhat[i][j] = newUs;
            }
        }
        setVelocityConstant(field);
        diff = sqrt(diff);
        DEBUGPRT("Solving U Momentum Eqn. Iter: %d, Diff: %lf\r\n", count, diff)
        count++;
    }

    // V
    diff = MAX;
    while((count<=MAXMEQNITER)&&(diff>CONV)){ // Iterate only inner nodes
        diff = 0;
        for(int i=0;i<field->N+2;i++){
            for(int j=0;j<field->M+1;j++){
                double a[6];
                momentum_equation_coeff(field,i,j,a,NOTHAT,ISV);
                double b = a[5]*get_v(field,i,j,NOTHAT);
                double newVs = (a[0]*get_v(field,i+1,j,HAT) 
                             + a[1]*get_v(field,i-1,j,HAT)
                             + a[2]*get_v(field,i,j+1,HAT)
                             + a[3]*get_v(field,i,j-1,HAT)
                             + b + field->dx*(get_P(field,i-1,j-1,field->P)-get_P(field,i-1,j,field->P))
                             + friction(field,i,j,ISV))/a[4];
                diff += pow(newVs-field->Vhat[i][j],2);
                field->Vhat[i][j] = newVs;
            }
        }
        setVelocityConstant(field);
        diff = sqrt(diff);
        DEBUGPRT("Solving V Momentum Eqn. Iter: %d, Diff: %lf\r\n", count, diff)
        count++;
    }

}

void SolvePressureCorrection(RFF* field, double** Pcorr){
    int count = 0;
    double diff = MAX;
    makezero(Pcorr,field->M,field->N);
    while((count<=MAXPCORITER)&&(diff>CONV)){
        diff = 0;
        for(int i=0;i<field->N;i++){
            for(int j=0;j<field->M;j++){
                double a[5];
                pressure_equation_coeff(field,i,j,a);
                double b = field->dy*(field->rho*get_u(field,i  ,j,HAT)
                                    - field->rho*get_u(field,i+1,j,HAT))
                         + field->dx*(field->rho*get_v(field,i+1,j,HAT)
                                    - field->rho*get_v(field,i+1,j+1,HAT));
                double newPs = (a[0]*get_P(field,i+1,j,Pcorr) 
                             + a[1]*get_P(field,i-1,j,Pcorr)
                             + a[2]*get_P(field,i,j+1,Pcorr)
                             + a[3]*get_P(field,i,j-1,Pcorr)
                             + b)/a[4];
                diff += pow(newPs-Pcorr[i][j],2);
                Pcorr[i][j] = newPs;
            }
        }
        diff = sqrt(diff);
        DEBUGPRT("Solving Pressure Eqn. Iter: %d, Diff: %lf\r\n",count,diff)
        count++;
    }
}

// Return the difference between Un and Un+1
double VelocityCorrection(RFF* field, double** Pcorr){
    double diff = 0;
    // U first
    for(int i=0;i<field->N+1;i++){
        for(int j=0;j<field->M;j++){
            double a[6];
            momentum_equation_coeff(field,i,j,a,NOTHAT,ISU);
            double newU = get_u(field,i,j,HAT) + (field->dy*(get_P(field,i-1,j,Pcorr)-get_P(field,i,j,Pcorr)))/a[4];
            diff += pow(newU-field->U[i][j],2);
            field->U[i][j] = newU;
        }
    }
    // V
    for(int i=0;i<field->N+2;i++){
        for(int j=0;j<field->M+1;j++){
            double a[6];
            momentum_equation_coeff(field,i,j,a,NOTHAT,ISV);
            double newV = get_v(field,i,j,HAT) + (field->dx*(get_P(field,i-1,j-1,Pcorr)-get_P(field,i-1,j,Pcorr)))/a[4];
            diff += pow(newV-field->V[i][j],2);
            field->V[i][j] = newV;
        }
    }
    setVelocityConstant(field);
    diff = sqrt(diff);
    return diff;
}
