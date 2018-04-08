#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>
#include <Windows.h>
#include <algorithm> 
#include <opencv2/opencv.hpp>
#include <sstream>
#include <stdio.h>

class Tools
{
public:
	Tools () {};
	~Tools () {};

	/*int z zakresu <min, max)*/
	static int r_int ( int min, int max );

	/*int z zakresu <0, 32767>*/
	static int r_int ();

	/*double z zakresu <min,max>*/
	static double r_db ( double min, double max );

	/*double z zakresu <0,1> */
	static double r_db ();

	bool init ();

	template <class type> static void print_tab ( type * tab, int size );

	template <class type> static void print_matrix ( type ** tab, int size_n, int size_m );

	template <class type> static std::stringstream tab_to_sstream ( type * tab, int size );

	template <class type> static std::stringstream matrix_to_sstream ( type ** tab, int size_n, int size_m );

	template <class type> static std::stringstream piec_d_to_sstream ( type ***** tab, int size_a, int size_b, int size_c, int size_d, int size_e );

	static void  piec_d_from_file (std::string adres, double ***** & tab, int& size_a, int &size_b, int &size_c, int &size_d, int &size_e );

	template <class type> static bool zapis_do_pliku ( std::string adres, type ob );

	template <class type> static bool odczyt_z_pliku ( std::string adres, type & ob );

	static std::vector <std::string> get_all_files_names_within_folder ( std::string folder );

	/*srednie rozwiazanie, duza zlozonosc obliczeniowa*/
	static std::vector <int> vektor_roznych_liczb ( int min, int max, int ilosc );

	/*zlozonosc obliczeniowa O( 3(max-min) ) ale zwraca tylko vektor o podanej dlugosci z zakresu*/
	static std::vector <int> vektor_roznych_liczb_ver2 ( int min, int max, int ilosc );

	//zlozonosc obliczeniowa O( 2(max-min) )
	static std::vector <int> vektor_roznych_liczb_ver2 ( int min, int max );

	static int * img_to_array ( std::string adres_katalogu, std::string nazwa_pliku, int & dlugosc );

	static int ** folder_of_img_to_array ( std::string adres_katalogu, int& ilosc_w_folderze, int & dlugosc_vektorow );

	static std::string int_to_fancy_vector ( int i );

	static int reverseInt(int i);

	static void read_mnist_data(std::string adres, int & ilosc_obrazow, int & wielkosc_obrazow, double ** & data);

	static void read_mnist_labels(std::string adres, int & ilosc_etykiet, double ** & etykiety);

	static void wagi_to_png(std::string nazwa, int ile_wag, int wys, int szer, double ** wagi);

	static void wagi_to_png_transp(std::string nazwa, int ile_wag, int wys, int szer, double ** wagi);
	

	/*Wymagany aby funkcje byly losowe*/
	static bool static_init ();
};

template<class type>
inline void Tools::print_tab ( type * tab, int size )
{
	for ( int i = 0; i < size; i++ )
	{
		std::cout << tab[i] << " ";
	}
	std::cout << "\n";
}

template<class type>
inline void Tools::print_matrix ( type ** tab, int size_n, int size_m )
{
	for ( int i = 0; i < size_n; i++ )
	{
		for ( int j = 0; j < size_m; j++ )
		{
			std::cout << tab[i][j] << " ";
		}
		std::cout << "\n";
	}
}

//pomianie 1 elementu na potrzeby latwiejszego dzialania
template<class type>
inline std::stringstream Tools::tab_to_sstream ( type * tab, int size )
{
	std::stringstream ss;
	//ss << size << "; ";
	for ( int i = 1; i < size; i++ )
	{
		ss << tab[i] << "; ";
	}
	ss << "\n";
	return ss;
}

template<class type>
inline std::stringstream Tools::matrix_to_sstream ( type ** tab, int size_n, int size_m )
{
	std::stringstream ss;
	ss << size_n << "; " << size_m << ";\n";
	for ( int i = 0; i < size_n; i++ )
	{
		for ( int j = 0; j < size_m; j++ )
		{
			ss << tab[i][j] << "; ";
		}
		ss << "\n";
	}
	return ss;
}

template<class type>
inline std::stringstream Tools::piec_d_to_sstream ( type ***** tab, int size_a, int size_b, int size_c, int size_d, int size_e )
{
	std::stringstream ss;
	ss << size_a << "; " << size_b << "; " << size_c << "; " << size_d << "; " << size_e << ";\n";
	for ( int i = 0; i < size_a; i++ )
	{
		for ( int j = 0; j < size_b; j++ )
		{
			for ( int k = 0; k < size_c; k++ )
			{
				for ( int l = 0; l < size_d; l++ )
				{
					for ( int m = 0; m < size_e; m++ )
					{
						ss << tab[i][j][k][l][m] << " ";
					}
				}
			}
		}
	}
	return ss;
}

inline void Tools::piec_d_from_file (std::string adres, double ***** & tab, int & size_a, int & size_b, int & size_c, int & size_d, int & size_e )
{
	std::ifstream ss;
	ss.open ( adres );
	ss >> size_a;
	ss.ignore ( 1 );
	ss >> size_b;
	ss.ignore ( 1 );
	ss >> size_c;
	ss.ignore ( 1 );
	ss >> size_d;
	ss.ignore ( 1 );
	ss >> size_e;
	ss.ignore ( 1 );
	tab = new double ****[size_a];
	for ( int i = 0; i < size_a; i++ )
	{
		tab[i] = new double ***[size_b];
		for ( int j = 0; j < size_b; j++ )
		{
			tab[i][j] = new double **[size_c];
			for ( int k = 0; k < size_c; k++ )
			{
				tab[i][j][k] = new double *[size_d];
				for ( int l = 0; l < size_d; l++ )
				{
					tab[i][j][k][l] = new double[size_e];
					for ( int m = 0; m < size_e; m++ )
					{
						ss >> tab[i][j][k][l][m];
						ss.ignore ( 1 );
					}
				}
			}
		}
	}
	ss.close ();
}



template<class type>
inline bool Tools::zapis_do_pliku ( std::string adres, type ob )
{
	std::ofstream plik;
	plik.open ( adres );
	if ( plik.is_open () )
	{
		plik << ob;
		plik.close ();
		std::cout << "Zapis do pliku: ";
		std::cout << adres;
		std::cout << '\n';
		return true;
	}
	else
	{
		std::cout << "Blad zapisu do: ";
		std::cout << adres;
		std::cout << '\n';
		return false;
	}
}

template<class type>
inline bool Tools::odczyt_z_pliku ( std::string adres, type & ob )
{
	std::ifstream plik;
	plik.open ( adres );
	if ( plik.is_open () )
	{
		plik >> ob;
		plik.close ();
		std::cout << "Wczytano z pliku:";
		std::cout << adres;
		std::cout << '\n';
		return true;
	}
	else
	{
		std::cout << "Blad wczytania z pliku: ";
		std::cout << adres;
		std::cout << '\n';
		return false;
	}
}

#endif // !TOOLS_H


