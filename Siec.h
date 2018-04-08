#ifndef SIEC_H
#define SIEC_H

#include "Tools.h"
#include <math.h>  
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>

class Siec
{
public:

	int n; //iloœæ wejsc
	int m; //ilosc neuronow warstwy wyjsciowej
	int l; //ilosc neuronow warstwy ukrytej

	double ** wagi_war_ukr;
	double ** wagi_war_wyj;

	double momentum;
	double ** zmiana_wagi_war_ukr;
	double ** zmiana_wagi_war_wyj;

	double * net_war_ukr;
	double * net_war_wyj;

	double * bias_war_ukr;
	double * bias_war_wyj;

	double * wyj_war_ukr;
	double * wyj_war_wyj;

	bool * neuron_aktywny;
	double prawd_wylaczenia=0;

	double * blad_war_ukr;
	double * blad_war_wyj;

	double * wejcia;
	double * oczekiwane_wyj;

	double war_uczenia;

	// ile wagi zostanie po redukcji
	double reg_L2= 0.999999;

	double ** ciagi_uczace;
	int ilosc_ciagow_uczacych = 0;
	double ** ciagi_uczace_war_oczekiwana;

	double ** ciagi_walidacyjne;
	int ilosc_ciagow_walidacyjnych = 0;
	double ** ciagi_walidacyjne_war_oczekiwana;

	double ** ciagi_testowe;
	int ilosc_ciagow_testowych = 0;
	double ** ciagi_testowe_war_oczekiwana;

	Siec ();
	~Siec ();

	/*N-Ilosc wejsc, M-ilosc neuronow wyjsciwych, L-ilosc neuronow warstwy ukrytej*/
	Siec ( int N, int M, int L, double War_ucz );
	bool reinit ( int N, int M, int L, double War_ucz );

	bool zmien_ilosc_neuronow_war_ukrytej ( int L );
	bool zmien_ilosc_neuronow_war_wyj ( int M );

	bool policz_net_ukr ();
	bool policz_net_wyj ();

	bool policz_wyj_ukr ();
	bool policz_wyj_wyj ();

	bool ustaw_wagi_war_ukr ( double min, double max );
	bool ustaw_wagi_war_wyj ( double min, double max );

	/*ustawia wszystkie wagi*/
	bool ustaw_wagi ( double min, double max );

	/*liczy wszystkie nety i wyj, nie ustawia wejsc*/
	double policz_wyj ();

	/*Przyblizenie = 0 wartosci zostana przyblizone do 1
	=1 przyblizona zostanie wartosc najwieksza do 1 reszta do 0
	=2 true wartosci*/
	double policz_wyjscie ( double * wejscie, int przyblizenie );

	double policz_blad_war_wyj ();
	double policz_blad_war_ukr ();

	double uaktualnij_wagi_wyj ();
	double uaktualnij_wagi_ukr ();

	double fun_binary_sigm ( double x );
	double fun_binary_sigm_p ( double x );
	double fun_ReLU ( double x );
	double fun_ReLU_p ( double x );

	/*0 - sigmoid, pp- ReLU*/
	int jaka_f_akt_ukr = 0;

	/*0- f akt softmax, pp- f akt sigmoid*/
	int jaka_f_akt_wyj = 0;

	double blad_sieci ();

	bool ustaw_wejscia ( double * wektor );
	bool ustaw_ocz_wyjsci ( double * wektor );

	double uczenie ( double * wektor_ucz, double * wektor_ocz );
	double epoki ( double ** wektor_ucz, double ** wektor_ocz, int ile_ciagow, int max_epok );

	double epoki_uczenie ( int max_epok, double zakres_bledu );

	/*Przechodzi jedna epoke i zwraca sume bledow uczenia*/
	double epoka ();

	bool wczytaj_ciag ( double ** & ciag_wejsciowy, int & ilosc, std::string adres );

	void pokaz_wagi ();

	double test ( double ** ciag_wejsciowy, double ** ciag_oczekiwany, int ilosc_ciagow );

	/*Blad srednio kwadratowy*/
	double test2 ( double ** ciag_wejsciowy, double ** ciag_oczekiwany, int ilosc_ciagow );

	/*przyblizenie = true: wynik 1 dla wartosci > prrogu,
	przyblizenie= false: wynik 1 dla wartosci najwiekszej
	zwraca true gdy ciagi sa "takie same"*/
	bool porownaj_wynik_z_ocz ( double * wynik, double * oczekiwany, int dlugosc, bool przyblizenie, double prog );

	bool save_to_csv ( std::string adres );

	/*Zwraca pare <srednia_ilosc_epok_nauczania, warjancja> */
	std::pair <double, double> badanie ( int ilosc_powotrzen, double wartosc_uczenia, double zakres_wag, double momentum, int jaka_funkcja_aktywacji, int liczba_neuronow_w_war_ukrytej );

	void badanie2 ( std::string adres,int max_epok, double wartosc_uczenia, double zakres_wag, double momentum, int jaka_funkcja_aktywacji, int liczba_neuronow_w_war_ukrytej );

	/*Hardcore */
	std::pair <double *****, double *****>  turbo_badanie ( int ilosc_powotrzen, std::vector<double> wartosc_uczenia, std::vector<double> zakres_wag, std::vector<double> momentum, std::vector<int> jaka_funkcja_aktywacji, std::vector<int> liczba_neuronow_w_war_ukrytej );

	void dawaj_skal_wagi_ukr(std::string adres);
	void dawaj_skal_wag_wyj(std::string adres);

	friend std::ostream & operator<< ( std::ostream& Strm_Wyj, const Siec & siec );
	friend std::istream & operator >> ( std::istream & Strm_We, Siec & siec );

};




#endif // !SIEC_H




