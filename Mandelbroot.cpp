# include <iostream>
# include <cstdlib>
# include <string>
# include <chrono>
# include <cmath>
# include <vector>
# include <fstream>
#include <mpi.h>

/** Une structure complexe est définie pour la bonne raison que la classe
 * complex proposée par g++ est très lente ! Le calcul est bien plus rapide
 * avec la petite structure donnée ci--dessous
 **/
struct Complex
{
    Complex() : real(0.), imag(0.)
    {}
    Complex(double r, double i) : real(r), imag(i)
    {}
    Complex operator + ( const Complex& z )
    {
        return Complex(real + z.real, imag + z.imag );
    }
    Complex operator * ( const Complex& z )
    {
        return Complex(real*z.real-imag*z.imag, real*z.imag+imag*z.real);
    }
    double sqNorm() { return real*real + imag*imag; }
    double real,imag;
};

std::ostream& operator << ( std::ostream& out, const Complex& c )
{
  out << "(" << c.real << "," << c.imag << ")" << std::endl;
  return out;
}

/** Pour un c complexe donné, calcul le nombre d'itérations de mandelbrot
 * nécessaires pour détecter une éventuelle divergence. Si la suite
 * converge, la fonction retourne la valeur maxIter
 **/
int iterMandelbrot( int maxIter, const Complex& c)
{
    Complex z{0.,0.};
    // On vérifie dans un premier temps si le complexe
    // n'appartient pas à une zone de convergence connue :
    // Appartenance aux disques  C0{(0,0),1/4} et C1{(-1,0),1/4}
    if ( c.real*c.real+c.imag*c.imag < 0.0625 )
        return maxIter;
    if ( (c.real+1)*(c.real+1)+c.imag*c.imag < 0.0625 )
        return maxIter;
    // Appartenance à la cardioïde {(1/4,0),1/2(1-cos(theta))}    
    if ((c.real > -0.75) && (c.real < 0.5) ) {
        Complex ct{c.real-0.25,c.imag};
        double ctnrm2 = sqrt(ct.sqNorm());
        if (ctnrm2 < 0.5*(1-ct.real/ctnrm2)) return maxIter;
    }
    int niter = 0;
    while ((z.sqNorm() < 4.) && (niter < maxIter))
    {
        z = z*z + c;
        ++niter;
    }
    return niter;
}

/**
 * On parcourt chaque pixel de l'espace image et on fait correspondre par
 * translation et homothétie une valeur complexe c qui servira pour
 * itérer sur la suite de Mandelbrot. Le nombre d'itérations renvoyé
 * servira pour construire l'image finale.
 
 Sortie : un vecteur de taille W*H avec pour chaque case un nombre d'étape de convergence de 0 à maxIter
 MODIFICATION DE LA FONCTION :
 j'ai supprimé le paramètre W étant donné que maintenant, cette fonction ne prendra plus que des lignes de taille W en argument.
 **/
void 
computeMandelbrotSetRow( int W, int H, int maxIter, int num_ligne, int* pixels)
{
    // Calcul le facteur d'échelle pour rester dans le disque de rayon 2
    // centré en (0,0)
    double scaleX = 3./(W-1);
    double scaleY = 2.25/(H-1.);
    //
    // On parcourt les pixels de l'espace image :
    for ( int j = 0; j < W; ++j ) {
       Complex c{-2.+j*scaleX,-1.125+ num_ligne*scaleY};
       pixels[j] = iterMandelbrot( maxIter, c );
    }
}

std::vector<int>
computeMandelbrotSet( int W, int H, int maxIter )
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::vector<int> pixels(W*H);
    start = std::chrono::system_clock::now();
    // On parcourt les pixels de l'espace image :
    for ( int i = 0; i < H; ++i ) {
      computeMandelbrotSetRow(W, H, maxIter, i, pixels.data() + W*(H-i-1) );
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Temps calcul ensemble mandelbrot : " << elapsed_seconds.count() 
              << std::endl;
    return pixels;
}


std::vector<int> computeMandparallele( int W, int H, int a, int b,int maxIter )//un efonction qui permet de devisier l'image ti chaque nombre de lignes va etre utilise par un processus
{
   
    std::vector<int> pixels(W*H);
   
    // On parcourt les pixels de l'espace image :
    for ( int i = a; i < b; ++i ) {
      computeMandelbrotSetRow(W, H, maxIter, i, pixels.data() + W*(i-a) );
    }
    
    return pixels;
}





/** Construit et sauvegarde l'image finale **/
void savePicture( const std::string& filename, int W, int H, const std::vector<int>& nbIters, int maxIter )
{
    double scaleCol = 1./maxIter;//16777216
    std::ofstream ofs( filename.c_str(), std::ios::out | std::ios::binary );
    ofs << "P6\n"
        << W << " " << H << "\n255\n";
    for ( int i = 0; i < W * H; ++i ) {
        double iter = scaleCol*nbIters[i];
        unsigned char r = (unsigned char)(256 - (unsigned (iter*256.) & 0xFF));
        unsigned char b = (unsigned char)(256 - (unsigned (iter*65536) & 0xFF));
        unsigned char g = (unsigned char)(256 - (unsigned( iter*16777216) & 0xFF));
        ofs << r << g << b;
    }
    ofs.close();
}



int main(int argc, char *argv[] ) 
 { 
    const int W = 800;
    const int H = 600;
    int numtasks , rank;
    //const int maxIter = 16777216;
    const int maxIter = 8*65536;
   MPI_Init ( &argc , &argv );//Initialisation de MPI
   MPI_Comm_size ( MPI_COMM_WORLD , & numtasks );//determine le nobre des taches
   MPI_Comm_rank ( MPI_COMM_WORLD , &rank);//lire le rang

   MPI_Status status ;


   //question 1
   int x=H/numtasks;
  
 /*std::chrono::time_point<std::chrono::system_clock> start, end;
 start = std::chrono::system_clock::now();
auto iters=computeMandparallele( W,  H, rank*x, (rank+1)*x, maxIter );
if(rank==0){
    //std::vector<int> iters=computeMandparallele( W,  H, rank*x, (rank+1)*x, maxIter );
//std::vector<int>result;
int i,j;
double scaleCol = 1./maxIter;//16777216
    std::ofstream ofs( "mandelbrot1.tga", std::ios::out | std::ios::binary );
    ofs << "P6\n"
        << W << " " << H << "\n255\n";*/



/*
   /** cettte patie n'est pas incluse juste pour un essai**/  

//std::vector<int> iters=computeMandparallele( W,  H, rank*x, (rank+1)*x, maxIter );
//iters.insert(result.end(), iters.begin(), iters.end());
//savePicture("mandelbrot1.tga", W, x, iters, maxIter);
   //for ( j = 1; j < numtasks ; j++){
/*for(j=0;j<iters.size();j++){
    result.push_back(iters[j]);
}*/

//MPI_Recv (&iters[0] , W * x, MPI_INT , k, 0, MPI_COMM_WORLD ,& status );
/*for(j=0;j<iters.size();j++){//j'ai essaye d'utiliser la fonction sans la modifier en utilant un vecteur sui s'ajoute a chaque fois mais pour le reste des proc s affiche en noir
    result.push_back(iters[j]);
}
}
savePicture("mandelbrot1.tga", W, H, result, maxIter);
*/
 
/***fin essai **/


/**suite**/
        /*for (  i = 0; i < W * x; ++i ) {
                double iter = scaleCol*iters[i];
                unsigned char r = (unsigned char)(256 - (unsigned (iter*256.) & 0xFF));
                unsigned char b = (unsigned char)(256 - (unsigned (iter*65536) & 0xFF));
                unsigned char g = (unsigned char)(256 - (unsigned( iter*16777216) & 0xFF));
                ofs << r << g << b;
            }
        for ( int k = 1; k <numtasks ; k++){
           
            MPI_Recv (&iters[0] , W *x, MPI_INT , k, 0, MPI_COMM_WORLD ,& status );
            
            for ( int i = 0; i < W * x; ++i ) {
                double iter = scaleCol*iters[i];
                unsigned char r = (unsigned char)(256 - (unsigned (iter*256.) & 0xFF));
                unsigned char b = (unsigned char)(256 - (unsigned (iter*65536) & 0xFF));
                unsigned char g = (unsigned char)(256 - (unsigned( iter*16777216) & 0xFF));
                ofs << r << g << b;
            }
        }
        ofs.close();
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Temps calcul ensemble mandelbrot : " << elapsed_seconds.count() << std::endl;
    }else{
        
        MPI_Send(&iters[0], W * x, MPI_INT, 0, 0, MPI_COMM_WORLD);
      
    }*/
/**fin question 1**/
    



//question 2 maitre esclave**/
std::chrono::time_point<std::chrono::system_clock> start, end;
 start = std::chrono::system_clock::now();

if ( rank == 0 ) // rank == 0 => master



    {
        
       
         std::vector<int>result;
        int count_task = 0;
        for ( int i = 1 ; i < numtasks; ++i ) {
            MPI_Send(&count_task, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            
            count_task += 1;
        }
        while ( count_task < H) {//numtasks==H
           //travail demandé

            std::vector<int> buff(W);
           
            MPI_Recv(&buff[0] , W , MPI_INT , MPI_ANY_SOURCE, 0, MPI_COMM_WORLD ,& status );
            for (int k = 0; k < W ; k++){
               
                result.push_back(buff[k]);
            }
          
            MPI_Send(&count_task, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
            
            count_task += 1 ;
        }
       
        count_task = -1;
        for ( int i = 1 ; i <numtasks ; ++i ){
            MPI_Send(&count_task, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
           
        }
        savePicture("mandelbrot.tga", W, H,  result, maxIter );
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Temps calcul ensemble mandelbrot : " << elapsed_seconds.count() << std::endl;





    }else{
    
        int num_task = 0;
      
        while (num_task != -1){
            MPI_Recv(&num_task, 1 , MPI_INT , MPI_ANY_SOURCE, 0, MPI_COMM_WORLD ,& status );
           
            if (num_task >= 0) {
             //travail demandé

                int num_ligne = num_task;//chaqsue ligne va etre traite par un proc jusqu'a la fin des lignes 
                std::vector<int> pi(W);
                computeMandelbrotSetRow(W, H, maxIter, num_ligne, pi.data());
                
               
                MPI_Send(&pi[0], W, MPI_INT, 0, 0, MPI_COMM_WORLD);
              
            }
        }

    }


MPI_Finalize();
    return EXIT_SUCCESS;
 }
    
