//tetrascatt_OK.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// private function
double ProductoMixto(  arma::vec a,  arma::vec b,  arma::vec c){
    double ret = a[ 0 ] * ( b[ 1 ] * c[ 2 ] - b[ 2 ] * c[ 1 ] );
           ret = ret + a[ 1 ] * ( b[ 2 ] * c[ 0 ] - b[ 0 ] * c[ 2 ] );
           ret = ret + a[ 2 ] * ( b[ 0 ] * c[ 1 ] - b[ 1 ] * c[ 0 ] );
  return(ret);
}

//Para adaptar poducto mixto a tetrascatt con armadillo, necesito que las cosas
// sean columnas.

// private function
double ProductoMixtoCol(  arma::mat a,  arma::mat b,  arma::mat c){
  double ret = a(0,0) * ( b( 1,0 ) * c(2,0) - b( 2,0 ) * c( 1,0 ) );
  ret = ret + a( 1,0 ) * ( b( 2,0 ) * c( 0,0 ) - b( 0,0 ) * c( 2,0 ) );
  ret = ret + a( 2,0 ) * ( b( 0,0 ) * c( 1,0 ) - b( 1,0 ) * c( 0,0 ) );
  return(ret);
}



// outdated function
arma::cx_vec tetrascatt_c_old(double cw, double g, double h,
                        arma::vec freq,arma::umat Tet, arma::mat Ver, arma::vec kversor){
  double Gamma0, f, Kmax,jac,ecs;
  double kab,kac,kad,kbc,kbd,kdc,ka,kb,kc,kd;
  int Ntet,  Nfreq; //Nver,
  arma::vec k;
  arma::cx_vec finf;
  arma::cx_double Suma, im,suma,finf_numer,suma22;
  arma::cx_double Dika,Dikad,Dikab,Dikac,Dikdc;
  arma::cx_double EDika,EDikc,EDikdkc;
  arma::mat VA,VB,VC,VD;
  ecs = 1.0e-14;
  im = arma::cx_double(0,1);

  // Medium parameters
  //Gamma0 = ( 1 - g * h * h ) / ( g * h * h ) - ( g - 1 ) / g ;
  //Ntet = size( Tet, 1 ) ;
  //Nver = size( Ver, 1 ) ;
  //Nfreq = size( freq, 1 ) ;
  //ts = zeros( Float64, Nfreq ) ;


  Gamma0 = (1 - g * h * h) / (g * h * h) - (g - 1) / g ;
  Ntet = Tet.n_rows;
  //Nver = Ver.n_rows;
  Nfreq = freq.n_elem;
  finf = finf.zeros(Nfreq);

  //  ts = zeros( Float64, Nfreq ) ;
  //  for q = 1 : Nfreq
  //    f = freq[ q ] ;
  //  Kmax = 2.0 * pi * f / cw ;
  //  k = Kmax / h * kversor ;
  //  Suma = Complex( 0 ) ;
  for ( int q=0; q<Nfreq; ++q){
    f = freq[q];
    Kmax = 2.0 * arma::datum::pi * f / cw ;
    k = Kmax / h * kversor ;
    Suma =0;
    //  for i = 1 : Ntet
    //# VC)rtices del tetraedro actual
    //    VA = Ver[ Tet[ i, 1 ], : ] ;
    //  VB = Ver[ Tet[ i, 2 ], : ] ;
    //  VC = Ver[ Tet[ i, 3 ], : ] ;
    //  VD = Ver[ Tet[ i, 4 ], : ] ;
    for ( int i=0; i<Ntet; ++i){
      //# VC)rtices del tetraedro actual
      // hay que restar uno porque todo empieza
      // en cero.
      VA = Ver.row(Tet(i, 0)-1).t();
      VB = Ver.row(Tet(i, 1)-1).t();
      VC = Ver.row(Tet(i, 2)-1).t();
      VD = Ver.row(Tet(i, 3)-1).t();
      jac = std::abs(ProductoMixtoCol(VB-VA, VC-VA, VD-VA));
      kab = dot(k, VB-VA);
      kac = dot(k, VC-VA);
      kad= dot(k, VD-VA);
      kbd = dot( k, VD - VB );
      kbc = dot( k, VC - VB );
      kdc = dot( k, VC - VD ) ;
      ka = dot( k, VA ) ;
      kb = dot( k, VB ) ;
      kc = dot( k, VC ) ;
      kd = dot( k, VD ) ;

      //# Productillos
      //      kab = dot3Fast( k, VB - VA ) ;
      //      kac = dot3Fast( k, VC - VA ) ;
      //      kad = dot3Fast( k, VD - VA ) ;
      //      kbc = dot3Fast( k, VC - VB ) ;
      //      kbd = dot3Fast( k, VD - VB ) ;
      //      kdc = dot3Fast( k, VC - VD ) ;
      //      ka = dot3Fast( k, VA ) ;
      //      kb = dot3Fast( k, VB ) ;
      //      kc = dot3Fast( k, VC ) ;
      //      kd = dot3Fast( k, VD ) ;
      // El gran IF
      if( std::abs(kab)<ecs ){
        if ( std::abs(kad ) < ecs ){
          if ( std::abs(kac ) < ecs ){
            suma =1.0 / 6 * exp( 2.0 * im * ka ) ;
            //suma=0 //sacar esto
          }
          else{ //expresion2
            Dikac = 2.0 * im * kac ;
            suma = 1.0 / pow(Dikac,3)  * exp( 2.0 * im * kc ) -
              1.0 / ( Dikac ) * exp( 2.0 * im * ka ) *
              ( 1.0 / 2.0 + 1.0 / Dikac + 1.0 / pow(Dikac,2) ) ;
            //suma=0 //sacar esto
          }
        }
        else{  //expresion3
          if (std::abs(kac ) < ecs){
            Dikad = 2.0 * im * kad ;
            Dika = 2.0 * im * ka ;
            Dikdc = 2.0 * im * kdc ;
            suma = 1.0 / ( Dikad * Dikdc ) * exp( 2.0 * im * kc ) -
              1.0 / ( Dikad * Dikdc ) * exp( 2.0 * im * kd ) -
              1.0 /  pow(Dikad,2)* exp( Dika ) -
              1.0 / ( 2.0 * Dikad ) * exp( Dika ) ;
            //suma=0 //sacar esto
          }
          else if  (std::abs( kac- kad) < ecs){ //EXPRESION 4
                 Dikad = 2.0 * im * kad ;
                 suma = 1.0 / ( pow(Dikad,2) ) * ( exp( 2.0 * im * ka ) *
                 ( 1.0 + 2.0 / Dikad ) +
                 exp( 2.0 * im * kd ) * ( 1.0 - 2.0 / Dikad ) ) ;
          }
          else{ //expresion 5
            Dikad = 2.0 * im * kad ;
            Dikac = 2.0 * im * kac ;
            Dikdc = 2.0* im * kdc ;
            suma = ( 1.0 / Dikad ) * ( exp( 2.0 * im * kc ) *
              ( 1.0 / ( Dikad * Dikdc) - 1.0 / ( Dikad * Dikac ) - 1.0 / pow(Dikac,2) ) -
              exp( 2.0 * im * kd ) / ( Dikad * Dikdc ) +
              exp( 2.0 * im * ka ) *
              ( 1.0 / ( Dikad * Dikac ) + 1.0 / Dikac + 1.0 / pow(Dikac,2) ) ) ;
            //suma=0 //sacar esto
          }
        }
      }
      else{
        if (std::abs(kad) < ecs){
          Dikab = 2.0 * im * kab ;
          EDika = exp( 2.0 * im * ka ) ;
          if (std::abs(kac ) < ecs){ // expresion 6
              suma = -1.0/(pow(Dikab,2))*(EDika)-  1.0/(pow(Dikab,3))*(EDika)-
              1.0/(2.0*(Dikab))*(EDika)+  1.0/(pow(Dikab,3))*(exp(2.0* im *kb));
            //suma=0 //sacar esto
          }
          else if  (std::abs(kab- kac)<ecs){ // EXPRESION 7
            suma = 1.0/(pow(Dikab,2))*(EDika)+
              2.0/(pow(Dikab,3))*(EDika)+
              1.0/(pow(Dikab,2))*(exp(2.0* im *kb))-
              2.0/(pow(Dikab,2))*(exp(2.0* im *kc));
          }
          else{ // expresion 8
              suma = im / ( 8.0 * kab ) * (
              exp( 2.0 * im * ka ) / kac * ( 1.0 / kab + 2.0 *im + 1.0 / kac ) -
                exp( 2.0 * im * kb ) / ( kab * kbc ) +
                exp( 2.0 * im * kc ) / kac * ( -1.0 / kab + kac /( kab * kbc ) -
                1.0 / kac ) ) ;

          }
        }
        else if  (std::abs(kad- kab)<ecs){   //#aca consideramos kad igual a kab
          Dikab = 2.0 * im * kab ;
          if (std::abs(kac ) < ecs){ // expresion 9
            EDika = exp( 2.0 * im * ka ) ;
            suma = 2.0/(pow(Dikab,3))*(EDika)+
              1.0/(pow(Dikab,2))*(exp(2.0* im *kb))-
              1.0/(pow(Dikab,3))*(exp(2.0* im *kb))+
              1.0/(pow(Dikab,2))*(EDika)-
              1.0/(pow(Dikab,3))*(exp(2.0* im *kd));
          }
          else if  (std::abs(kac -kab)< ecs){ //expresion 10
            suma = -1.0/(pow(Dikab,3))*(exp(2.0 * im *ka))+
              1.0/((2.0*Dikab))*(exp(2.0 * im *kb))+
              1.0/(pow(Dikab,3))*(exp(2.0 * im *kc))-
              1.0/((pow(Dikab,2)))*(exp(2.0 * im *kd));
          }
          else{//expresion 11
             suma = -1.0/(pow(Dikab,2)*(2.0 * im *kac))*(exp(2.0 * im *ka))-
              1.0/((Dikab)*(2.0 * im *kbc))*(exp(2.0 * im *kb))-
              1.0/((Dikab)*pow(2.0 * im *kbc,2))*(exp(2.0 * im *kb))+
              1.0/((Dikab)*pow(2.0 * im *kbc,2))*(exp(2.0 * im *kc))-
              1.0/(pow(Dikab,2)*(2.0 * im *kdc))*(exp(2.0 * im *kc))+
              1.0/(pow(Dikab,2)*(2.0 * im *kac))*(exp(2.0 * im *kc))+
              1.0/(pow(Dikab,2)*(2.0 * im *kdc))*(exp(2.0 * im *kd));
          }
        }
        else{ // aca consideramos kad no nulo
          if (std::abs(kac ) < ecs){ //expresion 12
            Dikab = 2.0 * im * kab ;
            Dikad = 2.0 * im * kad ;
            EDika = exp( 2.0 * im * ka ) ;
            suma = 1.0/((Dikab)*(Dikad)*(2.0 * im *kbd))*(exp(2.0* im *kd))-
              1.0/((Dikab)*(Dikad)*(2.0 * im *kbd))*(EDika)-
              1.0/((Dikab)*(Dikab)*(2.0 * im *kbd))*(exp(2.0* im *kb))+
              1.0/((Dikab)*(Dikab)*(2.0 * im *kbd))*(EDika)-
              1.0/((Dikab)*pow(Dikad,2))*(exp(2.0* im *kd))+
              1.0/((Dikab)*pow(Dikad,2))*(EDika)+
              1.0/((Dikab)*(Dikad))*(EDika);
          }
          else if  (std::abs(kac - kab)<ecs && std::abs(kad-kab)>ecs){ //EXPRESION 13
            Dikab = 2.0 * im * kab ;
            suma = 1.0/((Dikab)*(2.0 * im *kdc)*(2.0 * im *kbd))*(exp(2.0 * im *kc))-
              1.0/((Dikab)*(2.0 * im *kdc)*(2.0 * im *kbd))*(exp(2.0 * im *kd))-
              1.0/((Dikab)*(2.0 * im *kbd))*(exp(2.0 * im *kb))-
              1.0/((Dikab)*(2.0 * im *kdc)*(Dikad))*(exp(2.0 * im *kc))+
              1.0/((Dikab)*(2.0 * im *kdc)*(Dikad))*(exp(2.0 * im *kd))+
              1.0/((Dikad)*pow(Dikab,2))*(exp(2.0 * im *kc))-
              1.0/((Dikad)*pow(Dikab,2))*(exp(2.0 * im *ka));
            }
          else if  (std::abs(kac - kad)<ecs && std::abs(kab-kac)>ecs){ //EXPRESION 14
            Dikab = 2.0 * im * kab ;
            suma = 1.0/((Dikab)*(2.0 * im *kbd))*(exp(2.0 * im *kd))-
              1.0/((Dikab)*(2.0 * im *kbc)*(2.0 * im *kbd))*(exp(2.0 * im *kc))+
              1.0/((Dikab)*(2.0 * im *kbc)*(2.0 * im *kbd))*(exp(2.0 * im *kb))-
              1.0/((Dikab)*(Dikad))*(exp(2.0 * im *kd))+
              1.0/((Dikab)*pow(Dikad,2))*(exp(2.0 * im *kc))-
              1.0/((Dikab)*pow(Dikad,2))*(exp(2.0 * im *ka));
          }
          else{ //EXPRESION 15
            EDikc = exp( 2.0 * im * kc ) ;
            EDikdkc = exp( 2.0 * im *( kd-kc ) ) ;
            suma = im * EDikc / ( 8.0 * kab ) * (1.0 / kbd * ( ( 1.0 -
            EDikdkc ) / kdc - ( 1.0 - exp( 2.0*im*( kb-kc ) ) ) / kbc ) +
            1.0/kad * ( ( 1.0 - exp( 2.0*im*( ka-kc ) ) )/kac -
            ( 1.0 - EDikdkc ) / kdc ) ) ;
          }
        }
      }


      Suma =  Suma + jac *suma* Kmax*Kmax;
    }
    finf_numer = Suma * Gamma0/ ( 4.0 * arma::datum::pi );;
    finf[q] = finf_numer;
    //  ts[q] = 20.0 * abs(Suma)  ;
    // ts[q] = 20.0 *  abs(suma22)  ;
  }
  arma::cx_vec ret = finf;
  return(ret);
}

///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////EMPIEZA READAPTACION ELEGANTE DE EDMUNDO///////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


double SinC(double t){
  //chequear con edmundo los valores...
  double ret;
  if (std::abs(t)<1.0e-14)
    return 1;
  ret = sin(t)/t;
  return ret;
}

arma::cx_double  ReevaluacionSinc(double x){
  // taylor expansion of (sin(x)exp(ix) - 1)/(2ix)
  arma::cx_double im;
  im = arma::cx_double(0,1);
  arma::cx_double ret;
  double ecs = 1e-01;
  if ( (std::abs(x)<ecs)){
    ret = 0.5 +  1.0 * im * x/3.0 -pow(x,2)/6.0 - 1.0* im * pow(x,3)/15.0 + pow(x,4)/45.0;
  }
  else {
    ret = (SinC(x) * exp(im * x) - 1.0)/(2.0*im*x);
  }

   return ret;
}







// function OrdenarVertices( kversor, V1, V2, V3, V4 )
// # Orden: 1 2 3 4
//   EcS = 1e-8 ;
// Comb = [ [V1,V2,V3,V4], [V2,V3,V4,V1], [V3,V4,V1,V2], [V4,V1,V2,V3] ] ;
// for OrdenV in Comb
//   kab = dot3Fast( kversor, OrdenV[2] - OrdenV[1] ) ;
// kad = dot3Fast( kversor, OrdenV[4] - OrdenV[1]  ) ;
// # COMENTADO EN CODIGO ORIGINAL			kbc = dot3Fast( kversor, OrdenV[3] - OrdenV[2]  ) ;
// if abs( kab ) > EcS && abs( kad ) > EcS
//   return OrdenV ;
// end
//   end
//   println("ERROR en v?rtices. Posible tetraedro con problemas: ", "v1 :", V1 ) ;
// println("v2 :", V2 ) ;
// println("v3 :", V3 ) ;
// println("v4 :", V4 ) ;

// end


arma::mat OrdenarVertices(arma::vec kversor, arma::vec V1, arma::vec V2,arma::vec V3, arma::vec V4 )
{
  //arma::Row<int> OrdenarVertices(arma::vec kversor, arma::vec V1, arma::vec V2,arma::vec V3, arma::vec V4 ){
  arma::Mat<int> comb(24,4,arma::fill::zeros);
  arma::mat aux(3,4,arma::fill::randu);
  arma::mat ret(3,4,arma::fill::randu);
  arma::vec ab,ad;
  arma::Row<int> idx;
  double kab, kad, ecs;

  aux.col(0) = V1;
  aux.col(1) = V2;
  aux.col(2) = V3;
  aux.col(3) = V4;

  // inculyo todos los ordenes relativos fijando un vertice.
  comb(0,0) = 0;
  comb(0,1) = 1;
  comb(0,2) = 2;
  comb(0,3) = 3;
  //
  comb(1,0) = 0;
  comb(1,1) = 1;
  comb(1,2) = 3;
  comb(1,3) = 2;
  //
  comb(2,0) = 0;
  comb(2,1) = 2;
  comb(2,2) = 1;
  comb(2,3) = 3;
  //
  comb(3,0) = 0;
  comb(3,1) = 2;
  comb(3,2) = 3;
  comb(3,3) = 1;
  //
  comb(4,0) = 0;
  comb(4,1) = 3;
  comb(4,2) = 1;
  comb(4,3) = 2;
  //
  comb(5,0) = 0;
  comb(5,1) = 3;
  comb(5,2) = 2;
  comb(5,3) = 1;
  //

  // inculyo todos los ordenes relativos fijando un vertice.
  comb(6,0) = 1;
  comb(6,1) = 0;
  comb(6,2) = 2;
  comb(6,3) = 3;
  //
  comb(7,0) = 1;
  comb(7,1) = 0;
  comb(7,2) = 3;
  comb(7,3) = 2;
  //
  comb(8,0) = 1;
  comb(8,1) = 2;
  comb(8,2) = 0;
  comb(8,3) = 3;
  //
  comb(9,0) = 1;
  comb(9,1) = 2;
  comb(9,2) = 3;
  comb(9,3) = 0;
  //
  comb(10,0) = 1;
  comb(10,1) = 3;
  comb(10,2) = 0;
  comb(10,3) = 2;
  //
  comb(11,0) = 1;
  comb(11,1) = 3;
  comb(11,2) = 2;
  comb(11,3) = 0;
  //

  // inculyo todos los ordenes relativos fijando un vertice.
  comb(12,0) = 2;
  comb(12,1) = 0;
  comb(12,2) = 1;
  comb(12,3) = 3;
  //
  comb(13,0) = 2;
  comb(13,1) = 0;
  comb(13,2) = 3;
  comb(13,3) = 1;
  //
  comb(14,0) = 2;
  comb(14,1) = 1;
  comb(14,2) = 0;
  comb(14,3) = 3;
  //
  comb(15,0) = 2;
  comb(15,1) = 1;
  comb(15,2) = 3;
  comb(15,3) = 0;
  //
  comb(16,0) = 2;
  comb(16,1) = 3;
  comb(16,2) = 0;
  comb(16,3) = 1;
  //
  comb(17,0) = 2;
  comb(17,1) = 3;
  comb(17,2) = 1;
  comb(17,3) = 0;
  //

  // inculyo todos los ordenes relativos fijando un vertice.
  comb(18,0) = 3;
  comb(18,1) = 0;
  comb(18,2) = 1;
  comb(18,3) = 2;
  //
  comb(19,0) = 3;
  comb(19,1) = 0;
  comb(19,2) = 2;
  comb(19,3) = 1;
  //
  comb(20,0) = 3;
  comb(20,1) = 1;
  comb(20,2) = 0;
  comb(20,3) = 2;
  //
  comb(21,0) = 3;
  comb(21,1) = 1;
  comb(21,2) = 2;
  comb(21,3) = 0;
  //
  comb(22,0) = 3;
  comb(22,1) = 2;
  comb(22,2) = 0;
  comb(22,3) = 1;
  //
  comb(23,0) = 3;
  comb(23,1) = 2;
  comb(23,2) = 1;
  comb(23,3) = 0;
  //




  ecs = 1.0e-14;

  for ( int i=0; i<24; ++i){

    //    kab = dot3Fast( kversor, OrdenV[2] - OrdenV[1] ) ;
    ab = aux.col(comb(i,1))- aux.col(comb(i,0));
    kab = dot(kversor,ab) ;


    //    kad = dot3Fast( kversor, OrdenV[4] - OrdenV[1]  ) ;
    ad = aux.col(comb(i,3)) - aux.col(comb(i,0));
    kad = dot(kversor,ad);

    //      if abs( kab ) > EcS && abs( kad ) > EcS
    //        return OrdenV ;
    //      end
    //     end


    if ( (std::abs( kab ) > ecs) && (std::abs( kad )) > ecs){

      //idx = comb.row(i).t();
      //ret=aux.cols(comb.row(i));
      int  h=0;
      for (int j =0; j<4;++j){
        ret.col(h)=aux.col(comb(i,j));
        h=h+1;
      }

      //arma::Row<int> ret1 = comb.row(i);
      return(ret);
    }
  }
//agrego para que no tire warning
  return(ret);
}

//' @title tetrascatt_c
//'
//' @description Computes scattering from a volumetric mesh
//' efficientlty, it is an auxiliary function called by tetrascatt function.
//' @seealso \code{\link{tetrascatt}}
//' @param  cw sound speed in the water in m/s
//'
//' @param  g  density constrast value, i.e  g= rho1/rhow, where rho1 and rhow
//' are the density values of the stcatterer and the unbounded medium
//' respectively
//'
//' @param  h  density sound speed constrast value   i.e  h= c1/cw, where c1 is
//' the sound speed of the stcatterer
//'
//' @param  freq an array of frequencies, where the scattering is computed.
//'
//' @param  Ver  a matrix with the vertex of the tetrahedra, each vertex has to
//' have three coordinates.
//' @param  Tet a matrix containing the four index of each tetrahedron.
//' @param  kversor  three component vector that indicates the direction of the
//' incident plane wave.
//' @return  A complex number array which contains the backward
//' differential far-field scattering cross-section (f infinity)
//' values at each frequency.
//' @export
// [[Rcpp::export]]
arma::cx_vec tetrascatt_c(double cw, double g, double h,
                          arma::vec freq,arma::umat Tet,
                          arma::mat Ver,
                          arma::vec kversor){
  double Gamma0, f, Kmax,jac,ecs;
  double kab,kac,kad,kbc,kbd,kdc,ka,kb,kc,kd;
  int Ntet,  Nfreq; //Nver,
  arma::vec k;
  arma::cx_vec finf;
  arma::cx_double Suma, im,suma,finf_numer,suma22,ExpKdSincExpKdc;
  arma::cx_double Dika,Dikad,Dikab,Dikac,Dikdc;
  arma::cx_double EDika,EDikc,EDikdkc;
  arma::mat VA,VB,VC,VD;
  ecs = 1.0e-14;
  im = arma::cx_double(0,1);
  arma::mat return_ordenar_vertices(3,4,arma::fill::randu);
  // Medium parameters
  //Gamma0 = ( 1 - g * h * h ) / ( g * h * h ) - ( g - 1 ) / g ;
  //Ntet = size( Tet, 1 ) ;
  //Nver = size( Ver, 1 ) ;
  //Nfreq = size( freq, 1 ) ;
  //ts = zeros( Float64, Nfreq ) ;


  Gamma0 = (1 - g * h * h) / (g * h * h) - (g - 1) / g ;
  Ntet = Tet.n_rows;
  //Nver = Ver.n_rows;
  Nfreq = freq.n_elem;
  finf = finf.zeros(Nfreq);

  //  ts = zeros( Float64, Nfreq ) ;
  //  for q = 1 : Nfreq
  //    f = freq[ q ] ;
  //  Kmax = 2.0 * pi * f / cw ;
  //  k = Kmax / h * kversor ;
  //  Suma = Complex( 0 ) ;
  for ( int q=0; q<Nfreq; ++q){
      f = freq[q];
      Kmax = 2.0 * arma::datum::pi * f / cw ;
      k = Kmax / h * kversor ;
      Suma =0;
    for ( int i=0; i<Ntet; ++i){

      // V?rtices del tetraedro actual
      //VA, VB, VC, VD = OrdenarVertices( kversor, Ver[ Tet[ i, 1 ], : ],
      //                                        Ver[ Tet[ i, 2 ], : ], Ver[ Tet[ i, 3 ], : ],
      //                                                                  Ver[ Tet[ i, 4 ], : ] ) ;




      // hay que restar uno porque todo empieza
      // en cero.
      VA = Ver.row(Tet(i, 0)-1).t();
      VB = Ver.row(Tet(i, 1)-1).t();
      VC = Ver.row(Tet(i, 2)-1).t();
      VD = Ver.row(Tet(i, 3)-1).t();

      return_ordenar_vertices = OrdenarVertices(kversor,VA,VB,VC,VD);
      VA = return_ordenar_vertices.col(0);
      VB = return_ordenar_vertices.col(1);
      VC = return_ordenar_vertices.col(2);
      VD = return_ordenar_vertices.col(3);

      jac = std::abs(ProductoMixtoCol(VB-VA, VC-VA, VD-VA));
      kab = dot(k, VB-VA);
      kac = dot(k, VC-VA);
      kad= dot(k, VD-VA);
      kbd = dot( k, VD - VB );
      kbc = dot( k, VC - VB );
      kdc = dot( k, VC - VD ) ;
      ka = dot( k, VA ) ;
      kb = dot( k, VB ) ;
      kc = dot( k, VC ) ;
      kd = dot( k, VD ) ;

      // Integraci?n sobre el tetrahedro
      //    if abs( kbd ) < ecs
      if (std::abs( kbd ) < ecs){
        suma = 1.0 / ( 2.0 * im * kab ) *
          ( exp( 2.0*im*kd ) * ReevaluacionSinc( kdc ) +
          //( exp( 2.0*im*kd ) * ReevaluacionSinc( -1.0* kdc ) +
          ( exp( im*(ka+kc) ) * SinC(kac) - exp( im*(kd+kc) ) * SinC(kdc) ) / (2.0*im*kad) );
          //IndicadorExpresiones[11] += 1;}
          //suma=0.0;
      }
      else{
        ExpKdSincExpKdc = exp( im*(kd+kc) ) * SinC(kdc) ;
        suma = -1.0 / ( 4.0 * kab ) *
          ( ( ExpKdSincExpKdc - exp( im*(kb+kc) ) * SinC(kbc) ) / kbd  +
          ( exp( im*(ka+kc) ) * SinC(kac) - ExpKdSincExpKdc ) / kad ) ;
        //IndicadorExpresiones[15] += 1 ;
      }

      Suma =  Suma + jac *suma* Kmax*Kmax;

    }

    finf_numer = Suma * Gamma0/ ( 4.0 * arma::datum::pi );
    finf[q] = finf_numer;
    //  ts[q] = 20.0 * abs(Suma)  ;
    // ts[q] = 20.0 *  abs(suma22)  ;
  }
  arma::cx_vec ret = finf;
  return(ret);
}


