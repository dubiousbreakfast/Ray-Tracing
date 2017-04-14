#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include "ReadNSCL.hh"
#include <Eigen/Core>
#include "/home/caleb/CppNumericalSolvers/include/cppoptlib/meta.h"
#include "/home/caleb/CppNumericalSolvers/include/cppoptlib/problem.h"
#include "/home/caleb/CppNumericalSolvers/include/cppoptlib/boundedproblem.h"
#include "/home/caleb/CppNumericalSolvers/include/cppoptlib/solver/bfgssolver.h"
#include "/home/caleb/CppNumericalSolvers/include/cppoptlib/solver/neldermeadsolver.h"
#include "/home/caleb/CppNumericalSolvers/include/cppoptlib/solver/lbfgsbsolver.h"
#include <gsl/gsl_histogram.h>

using namespace std;
using namespace cppoptlib;
using namespace Eigen;



//Construct histograms
MatrixXd histogram(VectorXd& x, int rebin_factor=1){

  
  //commonly used variables
  int min = x.minCoeff();
  int max = x.maxCoeff();
  int size = (max - min + 1)/rebin_factor; //scale bins
  
  // define bins
  size_t Bins = (size);
  double range[size+1]; //range has to be one more than bin count
  //create range array
  for (int i=0; i<size+1; i++){
    range[i] = min+i;
  }
  
  //Set up the gsl histogram
  gsl_histogram* hist = gsl_histogram_alloc(Bins);
  assert( hist != NULL );
  gsl_histogram_set_ranges(hist, range, size+1);

  //put the data in
  for (int i=0; i<x.rows(); i++){
    double ele = x(i);
    gsl_histogram_increment(hist,ele);
  }

  //convert to channel vs count matrix
  MatrixXd result(size,2);
  for (int i=0; i<Bins; i++){
    result(i,0) = range[i];
    result(i,1) = gsl_histogram_get(hist,i);
  }

  return result;
}


// define a collection of back and front position vectors

struct both_pos{
  VectorXd front;
  VectorXd back;
};




//Change the pesky std::vectors into Eigen vectors
VectorXd to_matrix(vector<int>& input){
  int size = input.size();
  VectorXd new_vector(size);//initialize the Eigen vector
  for (int i=0; i<size;i++){
    new_vector(i) = (float)input.at(i); //fill vector with std::vector elements
  }
  return new_vector;
}

//overload for floating vectors
VectorXd to_matrix(vector<double>& input){
  int size = input.size();
  VectorXd new_vector(size);//initialize the Eigen vector
  for (int i=0; i<size;i++){
    new_vector(i) = input.at(i); //fill vector with std::vector elements
  }
  return new_vector;
}



//sorting
both_pos apply_gate(VectorXd& f, VectorXd& b, double start, double stop){

  both_pos the_values;
  vector<double> gf; //gated front
  vector<double> gb; //gated back
  //loop through and apply gate conditions
  for (int i=0; i<f.rows(); i++){
    //Only gate on front position though
    if (f(i)>start && f(i)<stop && b(i)!=0)
      {
	gf.push_back(f(i));
	gb.push_back(b(i));
      }
  }
  //finally put both into a both_pos struct
  the_values.front = to_matrix(gf);
  the_values.back = to_matrix(gb);
  return the_values;
}


//trying to chop of tails with a bunch of zeros

MatrixXd chop_tails(MatrixXd peak){
  

  double max_count = peak.col(1).maxCoeff();
  int chop_start = 1;
  int chop_stop = 1;
  
  //choose our stop value between 10% of peak max or 10 counts
  double cut_level;
  if (max_count*.01 > 5.0){
    cut_level = max_count*.01;
  }
  else{
    cut_level = 5.0;
  }

  //now loop until the values get sufficiently "peaky"
  for(int i=0; i<peak.rows();i++){
    if (peak(i,1) > cut_level){
      chop_start = i;
      break;
    }
  }
  //get the stop of the peak
  for(int i=(peak.rows()-1); i<peak.rows(); i--){
    if (peak(i,1) > cut_level){
      chop_stop = i;
      break;
    }
  }

  int size = (int)chop_stop - (int)chop_start; //number of rows needed
  MatrixXd chopped_peak = MatrixXd::Zero(size,2); //this will be our result
  
  //construct the new peak
  for(int i=0; i<chopped_peak.rows(); i++){
    int j = i + (int)chop_start; //index for original peak matrix
    chopped_peak(i,0) = peak(j,0);
    chopped_peak(i,1) = peak(j,1);
   }

  
  return chopped_peak;
  

}



//class for the Gaussian fits.

class chi_square : public Problem<double> {

private:
  VectorXd data_x;
  VectorXd data_y;
  double Area;
     
public:
  using typename cppoptlib::Problem<double>::Scalar;
  using typename cppoptlib::Problem<double>::TVector;
    

  double value(const TVector &x) {
    const TVector r1 = sqrt(Area/(pow(x[1],2)*2.0*M_PI))*exp(-.5*pow((data_x.array()-x[0]),2.0)/pow(x[1],2));
    const TVector r2 = pow((data_y.array() - r1.array()),2)/(data_y.array());
    return r2.sum();    
  }
  
  
  void set_data(VectorXd& x, VectorXd& y,double A){
    data_x = x; //channel number
    data_y = y; //counts in channel
    Area = A; //total number of counts
  }
  
};



//Fit Gaussian for a histogram and return the FWHM

double gaus_fit(MatrixXd& histo){
  VectorXd x = histo.col(0);
  VectorXd y = histo.col(1);
  chi_square f;
  f.set_data(x,y,y.sum());
  BfgsSolver<chi_square> solver;
  VectorXd mu_sigma(2);
  mu_sigma[0] = x.minCoeff()+(x.maxCoeff()-x.minCoeff())/2.0; //initial mu
  mu_sigma[1] = 50.0; //initial sigma^2
  solver.minimize(f,mu_sigma);
  return abs(mu_sigma[1])*2*sqrt(2*log(2));
}



//Find the reconstructed focal plane given the parameters 
VectorXd ray_trace(VectorXd P1, VectorXd P2, double H, double alpha){
  double S = 10.0; //parameter from detector in cm
  VectorXd num = (H*(P1.array()-P2.array()))/cos(alpha)+(S*P2.array());
  VectorXd dem = S+(tan(alpha)*(P1.array()-P2.array()));
  VectorXd pos_focal_plane = round(num.array()/dem.array());
  return pos_focal_plane;
}



//Find the reconstructed focal plane given the parameters 
VectorXd linear_ray_trace(VectorXd P1, VectorXd P2, double H){
  double S = 10.0; //parameter from detector in cm
  //taken from our paper 
  VectorXd pos_focal_plane = H/S*(P1-P2)+P2;
  return pos_focal_plane;
}



//class for the Focal Plane fits.

class FP_Fit : public Problem<double> {

private:
  vector<both_pos> peaks;

  
public:
  using typename cppoptlib::Problem<double>::Scalar;
  using typename cppoptlib::Problem<double>::TVector;

  
  
  
  void get_data(vector<both_pos>& all){
    peaks = all;
  }
  
  double value(const TVector& x){
    const double t1 = prelude_to_the_objective(x[0]);
    return t1;
  }

  double prelude_to_the_objective(double a){
    double total_sigma = 0.0;
    for (int i=0; i<peaks.size();i++){
      VectorXd at_plane = linear_ray_trace(peaks.at(i).front,peaks.at(i).back,a); //construct new focal plane
      int rebin = 1; //factor that we will increment to rebin histogram to remove bins with zero counts
      bool a_zero = true;
      MatrixXd new_peak = histogram(at_plane,rebin); //histogram it
      new_peak = chop_tails(new_peak); //chop of tails that might have zero counts before the rebin
      while (a_zero == true){
	bool zeros = (new_peak.col(1).array() == 0).any(); //check if there are zero count bins
	if (zeros){
	  rebin = rebin*2; //rebin by even numbers
	  new_peak = histogram(at_plane,rebin);
	}
	else{
	  a_zero = false; //check if there are empty bins break loop if there are none
	}
				    
      }

      double temp = gaus_fit(new_peak)*(double)rebin; //adjust sigma so that rebin is accounted for 
      total_sigma = total_sigma+temp;
      cout << a/10.0 << " " << total_sigma << " " << rebin << endl;
    }
    
    return total_sigma;
  }
  
};



//read our input file that defines region of interest on the front position section

MatrixXd ReadGateFile(const char* filename){
  int Nrows; //first line of file will tell us number of gates
  
  MatrixXd input(1,2);
  ifstream datafile;
  datafile.open(filename);

  if (!datafile.is_open()){
    cout << "ERROR: There was a problem opening the file!" << endl;
    return input;
  }

  else{
    int i=0;
    //loop through file
    for(string line; getline(datafile,line);i++){
      if (i==0){
	Nrows = atoi(line.c_str());
	input.resize(Nrows,2);
      }
      else{
	//separate line by white space and append to matrix 
        double a,b;
	stringstream s(line);
	s >> a >> b;
	if (i <= Nrows){
	  input(i-1,0) = a;
	  input(i-1,1) = b;
	}
      }
    }

  }
  return input;
}  



//If the eye is used with care and appropriate knowledge, it will see this is just ReadNSCL.cpp copied over.

void ReadEventFile(const char* filename, vector<Event> *allEvents, int nEvents){

  ifstream datafile;
  datafile.open (filename, ios::binary);
  //datafile.open ("/home/longland/midas/online/data/v1730/run00012.mid", ios::binary);

  if (!datafile.is_open()){
    cout << "ERROR: There was a problem opening the file!" << endl;
    return;
  }

  cout << "----------------------------------------" << endl;

   // Read the Begin-of-run event
  Event tEvent;
  tEvent.ReadEvent(datafile); 

  allEvents->reserve(1000);
  // Make an event and append the allEvents vector
  //  int i=0;
  if(nEvents > 0){
    for(int i=0; i<nEvents; i++){
      Event e;
      int ret = e.ReadEvent(datafile);
      if(ret == END_RUN)break;           // Return when End-of-Run is found
      //if(ret != PHYSICS_EVENT)continue;  // If an event isn't a physics event, just skip it
      allEvents->push_back(e);
    }
  } else {
    while(true){
      Event e;
      int ret = e.ReadEvent(datafile);
      //cout << "E " << i++ << " ret = " << ret << endl;
      if(ret == END_RUN)break;           // Return when End-of-Run is found
      //if(ret != PHYSICS_EVENT)continue;  // If an event isn't a physics event, just skip it
      allEvents->push_back(e);
    }
  }
  
  datafile.close();  // Close the input file
}


int Event::ReadEvent(ifstream& file){

  char hbuf[4];

  file.read( hbuf, sizeof(hbuf) );  
  evtHead.Size = (int)hbuf[0] & 0xFF;
  
  file.read( hbuf, sizeof(hbuf) );
  evtHead.Serial = (int)hbuf[0] & 0xFF;
  //PrintEventHeader();

  // Make a vector of Data
  uint32_t dataSize=evtHead.Size-8;

  if(evtHead.Serial == BEGIN_RUN){
    char word32[4];
    char word8[1];
    file.read( word32, sizeof(word32));
    cout << "Run Number: " << (uint32_t)word32[0] << endl;
    file.read( word32, sizeof(word32)); // time offset
    file.read( word32, sizeof(word32)); // time stamp
    uint32_t nRead = dataSize-12;
    bool doprint=true;
    for(uint32_t i=0;i<nRead;i++){
      file.read( word8, sizeof(word8) );
      if(word8[0]=='\0')doprint=false;
      if(doprint)cout << word8[0];
    }
    cout << endl;
    return BEGIN_RUN;
  }
  if(evtHead.Serial == END_RUN){
    return END_RUN;
  }
  
  // Make a bank
  Bank b;
  if(evtHead.Serial == PHYSICS_EVENT && dataSize == 142){
    strncpy(b.header.Name,"ADC1",sizeof(b.header.Name));
    uint32_t nRead = dataSize/2;
    
    // Read the raw data into a buffer
    //    char databuf[2];
    uint16_t databuf2=0;
    rawData.resize(nRead);
    for(uint32_t i=0; i<nRead; i++){
      //file.read(databuf, sizeof(databuf));
      //      cout << i << endl;
      file.read((char *)&databuf2, sizeof(databuf2));
      //cout << hex << setw(4) << databuf2 << dec << " ";
      uint32_t reading = databuf2 & 0xFFF;
      uint32_t flag = (databuf2 & 0xF000)>>12;
      //cout << dec << flag << " ";
      //if((i+1)%8 == 0)cout << endl;
      if(flag == 6)reading=0;
      if(flag == 5)reading=0xFFF;
      rawData[i] = reading;
    }
    //cout << endl;

    // reorganize the data because it's weird
    b.Data.resize((rawData.size()-7)/2);
    uint32_t index = 3;
    for(uint32_t i=0; i<b.Data.size(); i++){
      b.Data[rawData[index+1]]=rawData[index];
      index += 2;
    }
    /*  
    for(int i=0; i<16; i++){
      cout << setw(5) << b.Data[i];
    }
    cout << endl;
    */
  } else if(evtHead.Serial == INCREMENTAL_SCALERS) {
    strncpy(b.header.Name, "SCAL", sizeof(b.header.Name));

    //cout << "SCAL" << endl;
    //   uint32_t tstart,tend,tstamp,nscalers;
    char databuf[4];
    file.read(databuf, sizeof(databuf));
    //tstart = (uint32_t)databuf[0] & 0xFFFFFFFF;
    //cout << (int)databuf[0] << endl;
    file.read(databuf, sizeof(databuf));
    //tend = (uint32_t)databuf[0] & 0xFFFFFFFF;
    file.read(databuf, sizeof(databuf));
    //tstamp = (uint32_t)databuf[0] & 0xFFFFFFFF;
    file.read(databuf, sizeof(databuf));
    //nscalers = (uint32_t)databuf[0] & 0xFFFFFFFF;

    uint32_t nRead = dataSize/4 - 4;
    b.Data.resize(nRead);
    for(uint32_t i=0; i<nRead; i++){
      file.read(databuf, sizeof(databuf));
      b.Data[i] = (uint32_t)databuf[0] & 0xFF ;
      //if((i % 10) == 0)cout << endl;
      //cout << setw(6) << " " << b.Data[i] ;
    }
    //cout << endl;
    
  } else {
    strncpy(b.header.Name,"UNKN",sizeof(b.header.Name));
    uint32_t nRead = dataSize/2;

    // Read the raw data into a buffer
    char databuf[2];
    rawData.resize(nRead);
    for(uint32_t i=0; i<nRead; i++){
      file.read(databuf, sizeof(databuf));
      rawData[i] = (uint16_t)databuf[0] & 0xFF;
    }
  }

  // Put this mank into the Banks vector
  Banks.push_back(b);
  
  return evtHead.Serial;
}



void Event::PrintEventHeader(){

  cout << "--- Event Header ---" << endl;
  cout << "Size: " << evtHead.Size << endl;
  cout << "Type: " << evtHead.Serial << " : " << TranslateType() << endl;

}

string Event::TranslateType(){

  if(evtHead.Serial == 1)return "Begin of Run";
  if(evtHead.Serial == 2)return "End of Run";
  if(evtHead.Serial == 3)return "Pause of Run";
  if(evtHead.Serial == 4)return "Resume of Run";
  if(evtHead.Serial == 20)return "Scaler Data";
  if(evtHead.Serial == 30)return "Physics Event";

  return "UNKNOWN TYPE";

}

void Event::PrintEventData(){

  vector<uint16_t> x = Data;
  vector<uint16_t>::const_iterator i;
  for(i = x.begin(); i != x.end(); ++i){
    cout << (*i) << "  ";
  }
  cout << endl;
}





int main(int argc, char ** argv){

  //The man himself

     std::cout << R"(
___  ____ ____ _  _ _    ____ _   _    ____ _  _ _  _ ___    _  _ ___     ____ _  _ ___      _ ____ _  _ _ ___   /
|__] |__| |__/ |_/  |    |___  \_/     [__  |__| |  |  |     |  | |__]    |__| |\ | |  \     | |__| |\/| |  |   / 
|__] |  | |  \ | \_ |___ |___   |      ___] |  | |__|  |     |__| |       |  | | \| |__/    _| |  | |  | |  |  . 

     )" << '\n';


  
  // Quit if filename not given
  if(argc < 2)
  {
    cout << "Usage: BSUaJ <NSCL File> <Gate-File> " << endl;
    return -1;
  }

  //get the file names
  string filename = argv[1];
  string gate_filename = argv[2];
  
  int nEvents = 0;

  
  // Make a big vector of every event
  vector<Event> allEvents;
  vector<Event>::const_iterator itEvent;

  // Read the events
  ReadEventFile(filename.c_str(), &allEvents, nEvents);
  
  
  // **********************************************************************
  // Everything after this point is the sort routine

  // ---INITIALISERS-------------------------------------------------------

  // Define all of the data channels
  vector<int> vPos1;
  vector<int> vPos2;
  vector<int> vDE;
  vector<int> vE;

  // Define temporary variables that will be used on an event by event basis


  // Scalers
  int sScalers[32]={0};
  int nScalerEvents=0;

  // ---Fill Raw Data-------------------------------------------------------
  for(itEvent = allEvents.begin(); itEvent != allEvents.end(); ++itEvent){
    // Get the event
    Event tEvent = (*itEvent);

    // get the banks stored in the event
    vector<Bank> x = tEvent.getBanks();
    vector<Bank>::const_iterator itBank;

    // Loop through the banks to find what I want
    for(itBank = x.begin(); itBank != x.end(); ++itBank){
      Bank tBank = (*itBank);

      // The ADC Bank
      if(strncmp(tBank.header.Name,"ADC1",4) == 0){
	vE.push_back(tBank.Data[0]);
	vDE.push_back(tBank.Data[1]);
	vPos1.push_back(tBank.Data[2]);
	vPos2.push_back(tBank.Data[3]);
      }

      // The Scaler Bank
      if(strncmp(tBank.header.Name,"SCAL",4) == 0){
	nScalerEvents++;
	for(int i=0; i<32; i++){
	  sScalers[i] += tBank.Data[i];
	}
      }
    } // Bank iteration

  } // Event iteration
  cout << "Boomshakalaka, read " << allEvents.size() << " events" << endl;

  //Start of the actual computation
  
  
  //Read in the gates
  MatrixXd Gates = ReadGateFile(gate_filename.c_str());

  cout << "The following gate regions have been defined" << endl;
  cout << "Channel Range " << Gates << endl;

  //convert to matrices
  VectorXd POS1 = to_matrix(vPos1); 
  VectorXd POS2 = to_matrix(vPos2);

  
  //Apply the gates to the position vectors to define the peak regions
  vector<both_pos> gPOS;
  for (int i=0; i<Gates.rows(); i++){
    double start = Gates(i,0);
    double stop = Gates(i,1);
    gPOS.push_back(apply_gate(POS1,POS2,start,stop)); //vector of all peak channels front and back
  }

  //Do the initial fits
  cout << "Our Initial Values are" << endl;
  double init_total_sigma = 0.0;
  for (int i=0; i<gPOS.size();i++){
    both_pos peaks = gPOS.at(i);
    MatrixXd temp_hist = histogram(peaks.front);
    temp_hist = chop_tails(temp_hist);
    bool zero = true;
    int rebin = 1; 
    while (zero){
      if ((temp_hist.col(1).array() == 0).any()){
	rebin = rebin*2;
	temp_hist = histogram(peaks.front,rebin);
      }
      else{
	zero = false;
      }
    }
    double sigma = gaus_fit(temp_hist);
    init_total_sigma = init_total_sigma+sigma;
    cout << "Peak " << i+1 << " Standard Deviation " << abs(sigma) << endl;
  }
  cout << "Initial total sigma: " << init_total_sigma << endl;

  //now do the ray trace fit
  FP_Fit Result;
  Result.get_data(gPOS);
  
  // //initial values
  double H = 10.0;
  double alpha = 20.0;
  //go to radians
  alpha = alpha*(M_PI/180);
  VectorXd h_alpha(1);
  h_alpha[0] = H;
  //set bounds
  // VectorXd lower(2);
  // lower[0] = 0.0;
  // lower[1] = 1.0;
  // VectorXd upper(2);
  // upper[0] = 100.0;
  // upper[1] = 85.0;

  // Result.setLowerBound(lower);
  // Result.setUpperBound(upper);
  //choose solver
  NelderMeadSolver<FP_Fit> fp_solver;

  //minimize
  fp_solver.minimize(Result,h_alpha);

  cout << "The best H value is: " << h_alpha[0] << endl;
  cout << "H/S is: " << h_alpha[0]/10.0 << endl;
  
  //Now we take best values and return the histogram
  VectorXd Best = linear_ray_trace(POS1,POS2,h_alpha[0]);
  MatrixXd Best_Hist = histogram(Best);
  cout << "Total Sigma " << Result(h_alpha) << endl;

  //write the file

  ofstream myfile;
  myfile.open ("best.dat");
  for (int i=0; i<Best_Hist.rows(); i++){
    myfile << Best_Hist(i,0)  << " " <<  Best_Hist(i,1) << endl;
  }
  
  myfile.close();
  
   
}


