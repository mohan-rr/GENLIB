/*! \file fondateur.h
\brief Interface des fonctions de simulation calcul de probabilite

Interface de toutes les fonctions en rapport avec le gene fondateur

\author S�bastien Leclerc 
\contributor Jean-Fran�ois Lefebvre

*/

#ifndef GENFOND
#define GENFOND

#include <RcppCommon.h>
#include <unordered_map>

int simul(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
		int lSimul, double* pdRetConj,double* pdRetSimul,double* pdRetProp,double* probRecomb,double probSurvieHomo,int printprogress);

std::string simulhaplo(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int lNAncetre,
		int lSimul, double* probRecomb, std::unordered_map<int,haplotype*> *hapRef, std::string WD);

int getNumberRec(double* probRecomb, int sex);

double getRandomNumber(int exponential);

//int descendreHaplotypes(CIndSimul* Ordre_tmp, double probHap); //, /**/std::unordered_map<std::string, haplotype*>/* const std::unordered_map<std::string, haplotype*> &*/*hapRef);
//void makeRecomb( CIndSimul* Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, double posRecomb, int& cle );
void makeRecombM( CIndSimul* Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, double* posRecomb, int& cle );
void makeRecombF( CIndSimul* Ordre_tmp, std::unordered_map<int, haplotype*> *hapRef, double probHap, int nbRecomb, double* posRecomb, int& cle );

void recombine( haplotype* hapBegin, haplotype* hapEnd, haplotype* hapChild, int nbRecomb, double *posRecomb );

bool reconstruct(std::string WD, const std::string &hapfilename, const std::string &simufilename,const std::string &SNPposfilename,const int &BPsize);

bool ancestralseq(const std::string &fileName, std::unordered_map<float, std::string> &haploseqs);

std::vector<int> readSNPpos(const std::string &fileName);

int simulsingle(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
			 int lSimul, double* pdRetour,int printprogress);

int simulsingleFreq(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
				int lSimul, double* pdRetour,int printprogress);

SEXP simulsingleFct(int* SGenealogie, int * proposant, int lproposant, int* SplAncetre, int* SplAncEtatAll1, int* SplAncEtatAll2, int SlNAncetre,
				int SlSimul, SEXP SfctSousGrp, int Sprintprogress);

SEXP simulsingleProb(int* SGenealogie, int* SplProposant, int SlNProposant, int* SplAncetre,int SlNAncetres, int* SplAncEtat,SEXP mtProb,
				 int SlSimul, int Sprintprogress);

SEXP prob(int* Genealogie, int* plProposant, int* plProEtat,int lNProposant, int* plAncetre, int* plAncEtat, int lNAncetre,
	    double* pdRetConj,double* pdRetSimul,int printprogress,int onlyConj);

int CoefApparentement(int* Genealogie, int* plProposant, int lNProposant, int* plAncetre,	double* pdRetour,int DuppDetection, int printprogress);

#endif



