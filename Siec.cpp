#include "Siec.h"



Siec::Siec ()
{
}


Siec::~Siec ()
{
}

/*
N-iloœæ wejsc
M-ilosc neuronow warstwy wyjsciowej
L-ilosc neuronow warstwy ukrytej
War_ucz - wartosc parametru uczenia
*/
Siec::Siec ( int N, int M, int L, double War_ucz )
{
	war_uczenia = War_ucz;
	n = N;
	m = M;
	l = L;

	wagi_war_ukr = new double*[l];
	wagi_war_wyj = new double*[m];

	momentum = 0;
	zmiana_wagi_war_ukr = new double*[l];
	zmiana_wagi_war_wyj = new double*[m];

	net_war_ukr = new double[l] ();
	net_war_wyj = new double[m] ();

	bias_war_ukr = new double[l] ();
	bias_war_wyj = new double[m] ();

	wyj_war_ukr = new double[l] ();
	wyj_war_wyj = new double[m] ();

	neuron_aktywny = new bool[l] ();

	blad_war_ukr = new double[l] ();
	blad_war_wyj = new double[m] ();

	wejcia = new double[N] ();
	oczekiwane_wyj = new double[m] ();

	for ( int i = 0; i < l; i++ )
	{
		wagi_war_ukr[i] = new double[n] ();
		zmiana_wagi_war_ukr[i] = new double[n] ();
	}
	for ( int i = 0; i < m; i++ )
	{
		wagi_war_wyj[i] = new double[l] ();
		zmiana_wagi_war_wyj[i] = new double[l] ();
	}

	Tools::static_init ();
}

bool Siec::reinit ( int N, int M, int L, double War_ucz )
{
	//TO  DO ?? ale moze byc nie potrzebne
	return false;
}

bool Siec::zmien_ilosc_neuronow_war_ukrytej ( int L )
{
	if ( L == l ) return false;

	l = L;
	delete[] bias_war_ukr;
	delete[] wyj_war_ukr;
	delete[] net_war_ukr;
	delete[] blad_war_ukr;
	delete[] zmiana_wagi_war_ukr;
	delete[] zmiana_wagi_war_wyj;
	delete[] wagi_war_ukr;
	delete[] wagi_war_wyj;
	delete neuron_aktywny;

	bias_war_ukr = new double[l] ();
	wyj_war_ukr = new double[l] ();
	net_war_ukr = new double[l] ();
	blad_war_ukr = new double[l] ();

	wagi_war_ukr = new double*[l];
	wagi_war_wyj = new double*[m];
	zmiana_wagi_war_ukr = new double*[l];
	zmiana_wagi_war_wyj = new double*[m];
	neuron_aktywny = new bool[l];

	for ( int i = 0; i < l; i++ )
	{
		wagi_war_ukr[i] = new double[n] ();
		zmiana_wagi_war_ukr[i] = new double[n] ();
	}

	for ( int i = 0; i < m; i++ )
	{
		wagi_war_wyj[i] = new double[l] ();
		zmiana_wagi_war_wyj[i] = new double[l] ();
	}
	return true;
}

bool Siec::zmien_ilosc_neuronow_war_wyj ( int M )
{
	if ( M == m ) return false;

	m = M;

	delete[] bias_war_wyj;
	delete[] wyj_war_wyj;
	delete[] net_war_wyj;
	delete[] blad_war_wyj;
	delete[] zmiana_wagi_war_wyj;
	delete[] wagi_war_wyj;
	delete[] oczekiwane_wyj;

	wagi_war_wyj = new double*[m];
	zmiana_wagi_war_wyj = new double*[m];
	net_war_wyj = new double[m] ();
	bias_war_wyj = new double[m] ();
	wyj_war_wyj = new double[m] ();
	blad_war_wyj = new double[m] ();
	oczekiwane_wyj = new double[m] ();

	for ( int i = 0; i < m; i++ )
	{
		wagi_war_wyj[i] = new double[l] ();
		zmiana_wagi_war_wyj[i] = new double[l] ();
	}
	return true;
}

bool Siec::policz_net_ukr ()
{
	for ( int i = 0; i < l; i++ )
	{
		net_war_ukr[i] = 0;
		for ( int j = 0; j < n; j++ )
		{
			net_war_ukr[i] += wagi_war_ukr[i][j] * wejcia[j];
		}
		net_war_ukr[i] += bias_war_ukr[i];
	}
	return true;
}

bool Siec::policz_net_wyj ()
{
	for ( int i = 0; i < m; i++ )
	{
		net_war_wyj[i] = 0;
		for ( int j = 0; j < l; j++ )
		{
			net_war_wyj[i] += wagi_war_wyj[i][j] * wyj_war_ukr[j];
		}
		net_war_wyj[i] += bias_war_wyj[i];
	}
	return true;
}

bool Siec::policz_wyj_ukr ()
{
	for ( int i = 0; i < l; i++ )
	{
		if ( Tools::r_db () >= prawd_wylaczenia )
		{
			if ( jaka_f_akt_ukr == 0 )
				wyj_war_ukr[i] = fun_binary_sigm ( net_war_ukr[i] );
			else
				wyj_war_ukr[i] = fun_ReLU ( net_war_ukr[i] );
			neuron_aktywny[i] = true;
		}
		else
		{
			wyj_war_ukr[i] = 0;
			neuron_aktywny[i] = false;
		}
	}
	return true;
}

bool Siec::policz_wyj_wyj ()
{
	if ( jaka_f_akt_wyj == 0 )
	{
		double net_max = *std::max_element ( net_war_wyj, net_war_wyj + m );
		double sum = 0;
		for ( int i = 0; i < m; i++ )
		{
			wyj_war_wyj[i] = exp ( net_war_wyj[i] - net_max );
			sum += wyj_war_wyj[i];
		}
		for ( int i = 0; i < m; i++ )
		{
			wyj_war_wyj[i] /= sum;
		}
	}
	else
	{
		for ( int i = 0; i < m; i++ )
		{
			wyj_war_wyj[i] = fun_binary_sigm ( net_war_wyj[i] );
		}
	}

	return true;
}

bool Siec::ustaw_wagi_war_ukr ( double min, double max )
{
	for ( int i = 0; i < l; i++ )
	{
		for ( int j = 0; j < n; j++ )
		{
			wagi_war_ukr[i][j] = Tools::r_db ( min, max );
			zmiana_wagi_war_ukr[i][j] = 0;
		}
		bias_war_ukr[i] = Tools::r_db ( min, max );
	}
	return true;
}

bool Siec::ustaw_wagi_war_wyj ( double min, double max )
{
	for ( int i = 0; i < m; i++ )
	{
		for ( int j = 0; j < l; j++ )
		{
			wagi_war_wyj[i][j] = Tools::r_db ( min, max );
			zmiana_wagi_war_wyj[i][j] = 0;
		}
		bias_war_wyj[i] = Tools::r_db ( min, max );
	}
	return true;
}

bool Siec::ustaw_wagi ( double min, double max )
{
	return ustaw_wagi_war_ukr ( min, max ) && ustaw_wagi_war_wyj ( min, max );
}

double Siec::policz_wyj ()
{
	policz_net_ukr ();
	policz_wyj_ukr ();
	policz_net_wyj ();
	policz_wyj_wyj ();

	return 0.0;
}

double Siec::policz_wyjscie ( double * wejscie, int przyblizenie )
{
	ustaw_wejscia ( wejscie );
	policz_wyj ();
	double wynik = -1.0;
	if ( przyblizenie == 0 ) //przyblizanie wartosci
	{
		for ( int i = 0; i < m; i++ )
		{
			if ( wyj_war_wyj[i] > 0.5 )
			{
				wynik = i;
				std::cout << 1 << " ";
			}
			else
				std::cout << 0 << " ";
		}
	}
	else if ( przyblizenie == 1 ) //prawidziwe wartosci
	{
		double max = wyj_war_wyj[0];
		wynik = 0;
		//std::cout << wyj_war_wyj[0] << " ";
		for ( int i = 1; i < m; i++ )
		{
			//std::cout << wyj_war_wyj[i] << " ";
			if ( max < wyj_war_wyj[i] )
			{
				max = wyj_war_wyj[i];
				wynik = i;
			}
		}
		for ( int i = 0; i < m; i++ )
		{
			if ( i == wynik )
				std::cout << 1 << " ";
			else
				std::cout << 0 << " ";
		}
	}
	else
	{
		for ( int i = 0; i < m; i++ )
		{
			std::cout << wyj_war_wyj[i] << " ";
		}
	}

	std::cout << '\n';
	return wynik;
}

double Siec::policz_blad_war_ukr ()
{
	for ( int i = 0; i < l; i++ )
	{
		double suma = 0;
		for ( int j = 0; j < m; j++ )
		{
			suma += blad_war_wyj[j] * wagi_war_wyj[j][i];
		}
		if ( jaka_f_akt_ukr == 0 )
			blad_war_ukr[i] = fun_binary_sigm_p ( net_war_ukr[i] )*suma;
		else
			blad_war_ukr[i] = fun_ReLU_p ( net_war_ukr[i] )*suma;
	}
	return 0.0;
}

double Siec::policz_blad_war_wyj ()
{
	for ( int i = 0; i < m; i++ )
	{
		if ( jaka_f_akt_wyj == 0 )
		{
			double sum = 0;
			for ( int j = 0; j < l; j++ )
			{
				sum += wyj_war_ukr[j] / l;
			}
			blad_war_wyj[i] = sum*( oczekiwane_wyj[i] - wyj_war_wyj[i] );
		}
		else
		{
			blad_war_wyj[i] = ( oczekiwane_wyj[i] - wyj_war_wyj[i] )*fun_binary_sigm_p ( net_war_wyj[i] );
		}

	}
	return 0.0;
}

double Siec::uaktualnij_wagi_wyj ()
{
	double korekcja;
	for ( int i = 0; i < m; i++ )
	{
		//czy tutaj gdy wylaczamy neuron jego waga ma sie nie zmieniac przez regularyzacje czy jednak tak? chyba nie powwina.
		for ( int j = 0; j < l; j++ )
		{
			korekcja = war_uczenia * blad_war_wyj[i] * wyj_war_ukr[j];
			wagi_war_wyj[i][j] = reg_L2* wagi_war_wyj[i][j];
			wagi_war_wyj[i][j] += korekcja;// +momentum*zmiana_wagi_war_wyj[i][j];
			zmiana_wagi_war_wyj[i][j] = korekcja;
		}

		bias_war_wyj[i] += war_uczenia*blad_war_wyj[i];
	}
	return 0.0;
}

double Siec::uaktualnij_wagi_ukr ()
{
	double korekcja;
	for ( int i = 0; i < l; i++ )
	{
		if ( neuron_aktywny[i] )
		{
			for ( int j = 0; j < n; j++ )
			{
				korekcja = war_uczenia*blad_war_ukr[i] * wejcia[j];
				wagi_war_ukr[i][j] = reg_L2*wagi_war_ukr[i][j];
				wagi_war_ukr[i][j] += korekcja;// +momentum*zmiana_wagi_war_ukr[i][j];
				zmiana_wagi_war_ukr[i][j] = korekcja;
			}
			bias_war_ukr[i] += war_uczenia*blad_war_ukr[i];
		}
	}
	return 0.0;
}


double Siec::fun_binary_sigm ( double x )
{
	return 1 / ( 1 + exp ( -x ) );
}

double Siec::fun_binary_sigm_p ( double x )
{
	return fun_binary_sigm ( x )*( 1 - fun_binary_sigm ( x ) );
}

double Siec::fun_ReLU ( double x )
{
	/*Dlaczego przy zastosoaniu ReLU to nie dziala?*/
	return x < 0 ? 0 : x;
}

double Siec::fun_ReLU_p ( double x )
{
	/*Dlaczego przy zastosoaniu ReLU to nie dziala?*/
	return x < 0 ? 0 : 1;
}

double Siec::blad_sieci ()
{
	double blad = 0;
	for ( int i = 0; i < m; i++ )
	{
		blad += ( oczekiwane_wyj[i] - wyj_war_wyj[i] ) * ( oczekiwane_wyj[i] - wyj_war_wyj[i] );
	}
	return ( 0.5 *blad );
}

bool Siec::ustaw_wejscia ( double * wektor )
{
	//std::cout << "Ustawianie wejsc\n";
	for ( int i = 0; i < n; i++ )
	{
		wejcia[i] = wektor[i];
		//std::cout << wejcia[i] << " ";
	}
	//std::cout << "\n";
	return true;
}

bool Siec::ustaw_ocz_wyjsci ( double * wektor )
{
	//std::cout << "Ustawianie wyjsc\n";
	for ( int i = 0; i < m; i++ )
	{
		oczekiwane_wyj[i] = wektor[i];
		//std::cout << oczekiwane_wyj[i] << " ";
	}
	//std::cout << "\n";
	return true;
}

double Siec::uczenie ( double * wektor_ucz, double * wektor_ocz )
{
	ustaw_wejscia ( wektor_ucz );
	ustaw_ocz_wyjsci ( wektor_ocz );

	policz_net_ukr ();
	policz_wyj_ukr ();
	policz_net_wyj ();
	policz_wyj_wyj ();
	policz_blad_war_wyj ();
	policz_blad_war_ukr ();
	uaktualnij_wagi_wyj ();
	uaktualnij_wagi_ukr ();

	return blad_sieci ();
}

double Siec::epoki ( double ** wektor_ucz, double ** wektor_ocz, int ile_ciagow, int max_epok )
{
	double epoka = 0;
	bool stop = false;
	int *ciag_uzyty = new int[ile_ciagow] ();
	int ciag;
	while ( epoka < max_epok )
	{
		epoka++;
		std::cout << "epoka: " << epoka << "\n";

		//zerowanie wektora uzytych ciagow uczacych
		for ( int i = 0; i < ile_ciagow; i++ )
		{
			ciag_uzyty[i] = 0;
		}

		stop = true;
		for ( int i = 0; i < ile_ciagow; i++ )
		{
			//wybieranie nie uzytego ciagu uczacego, tzn losowa ich kolejnosc wynikowo
			ciag = Tools::r_int ( 0, ile_ciagow );
			while ( i < ile_ciagow && ciag_uzyty[ciag] == 1 )
			{
				ciag = Tools::r_int ( 0, ile_ciagow );
			}
			ciag_uzyty[ciag] = 1;
			std::cout << ciag << " ";

			double blad = uczenie ( wektor_ucz[ciag], wektor_ocz[ciag] );

			//if ( blad > zakres_bledu || blad < -zakres_bledu )
			{
				stop = false;
			}
		}
		std::cout << "\n";
	}

	return epoka;
}

double Siec::epoki_uczenie ( int max_epok, double zakres_bledu )
{
	std::ofstream MyExcelFile;
	MyExcelFile.open ( "Wyniki/test_epok2.csv" );
	MyExcelFile << "epoka;suma bledow;blad srednio kwad walidacji;blad wal;blad ucz;\n";

	int ktora_epoka = 0;
	bool stop = false;
	double suma_bledow = 0;
	std::cout << std::setprecision ( 4 );
	std::cout << std::fixed;
	std::cout.width ( 6 );
	while ( ktora_epoka < max_epok && !stop )
	{
		//stop = true;
		ktora_epoka++;
		suma_bledow = epoka ();

		if ( suma_bledow > zakres_bledu )
		{
			//stop = false;
		}
		else
		{
			//stop = false;
			//std::cout << "suma_bledow < zakres_bledu\n";
		}
		double blad_walid = test ( ciagi_walidacyjne, ciagi_walidacyjne_war_oczekiwana, ilosc_ciagow_walidacyjnych );
		double blad_uczenia = test ( ciagi_uczace, ciagi_uczace_war_oczekiwana, ilosc_ciagow_uczacych );
		double blad_sr_walidacji = test2 ( ciagi_walidacyjne, ciagi_walidacyjne_war_oczekiwana, ilosc_ciagow_walidacyjnych );
		//if ( blad_uczenia > 3.0 )
		{
			//stop = false;
			//Tools::zapis_do_pliku ( "siec.txt", *this );
			//system ( "pause" );
			//stop = true;
		}
		//if ( blad_walid <= 8.0 )
		{
			//Tools::zapis_do_pliku ( "siec2.txt", *this );
			//std::cout << "niski blad walid: " << blad_walid << "\n";
			//system ( "pause" );
		}

		std::cout << "epoka: " << ktora_epoka << "\tsuma bledow prop: " << suma_bledow << "\tblad sr wal: " << blad_sr_walidacji
			<< "\tblad wal: " << blad_walid << "\tblad ucz: " << blad_uczenia << "\n";
		Tools::zapis_do_pliku ( "Wyniki/S/siec" + std::to_string ( ktora_epoka ) + ".txt", *this );
		//save_to_csv ( "Wyniki/S/siec" + std::to_string ( ktora_epoka ) + ".csv" );

		MyExcelFile << ktora_epoka << ";" << suma_bledow << ";" << blad_sr_walidacji << ";" << blad_walid << ";" << blad_uczenia << ";\n";

		if (blad_uczenia < 0.00166667)
			stop = true;
	}
	MyExcelFile.close ();
	return ktora_epoka;
}

double Siec::epoka ()
{
	std::vector <int> wektor = Tools::vektor_roznych_liczb_ver2 ( 0, ilosc_ciagow_uczacych );
	double suma_bledow = 0;
	int ciag = 0;
	//std::cout << "\n";
	for ( int i = 0; i < ilosc_ciagow_uczacych; i++ )
	{
		//std::cout << i << " ";
		ciag = wektor[i];
		suma_bledow += uczenie ( ciagi_uczace[ciag], ciagi_uczace_war_oczekiwana[ciag] ) / ilosc_ciagow_uczacych;
	}
	//std::cout << "\n";
	wektor.clear ();
	return suma_bledow;
}

bool Siec::wczytaj_ciag ( double ** & ciag_wejsciowy, int & ilosc, std::string adres )
{
	std::ifstream plik;
	plik.open ( adres );
	int size = 0;
	if ( plik.is_open () )
	{
		plik >> ilosc;
		plik.ignore ( 1 );
		plik >> size;
		plik.ignore ( 1 );
		//std::cout << ilosc << " " << size << '\n';
		ciag_wejsciowy = new double *[ilosc];
		for ( int i = 0; i < ilosc; i++ )
		{
			ciag_wejsciowy[i] = new double[size] ();
			for ( int j = 0; j < size; j++ )
			{
				plik >> ciag_wejsciowy[i][j];
				plik.ignore ( 1 );
				//std::cout<< ciag_wejsciowy[i][j]<<" ";
			}
			//std::cout << "\n";
		}
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


void Siec::pokaz_wagi ()
{
	std::cout << std::fixed << std::showpos << std::setprecision ( 2 );
	std::cout << *this;
	std::cout << std::noshowpos << std::defaultfloat << std::setprecision ( 0 );
}

double Siec::test ( double ** ciag_wejsciowy, double ** ciag_oczekiwany, int ilosc_ciagow )
{
	double blad_walid = 0;
	for ( int i = 0; i < ilosc_ciagow; i++ )
	{
		ustaw_wejscia ( ciag_wejsciowy[i] );
		policz_wyj ();
		if ( !porownaj_wynik_z_ocz ( wyj_war_wyj, ciag_oczekiwany[i], m, true, 0.5 ) )
		{
			blad_walid++;
		}
	}
	//std::cout << "blad " << blad_walid / (double) ilosc_ciagow * 100 << "%\n";
	return ( blad_walid / (double) ilosc_ciagow * 100 );
}

double Siec::test2 ( double ** ciag_wejsciowy, double ** ciag_oczekiwany, int ilosc_ciagow )
{
	double blad_walid = 0;
	for ( int i = 0; i < ilosc_ciagow; i++ )
	{
		ustaw_wejscia ( ciag_wejsciowy[i] );
		policz_wyj ();
		for ( int j = 0; j < m; j++ )
			blad_walid += 0.5*( ciag_oczekiwany[i][j] - wyj_war_wyj[j] )*( ciag_oczekiwany[i][j] - wyj_war_wyj[j] ) / ilosc_ciagow;
	}
	return blad_walid;
}

bool Siec::porownaj_wynik_z_ocz ( double * wynik, double * oczekiwany, int dlugosc, bool przyblizenie, double prog )
{
	if ( przyblizenie )
	{
		double max = wynik[0];
		int index = 0;
		for ( int i = 1; i < dlugosc; i++ )
		{
			if ( max < wynik[i] )
			{
				max = wynik[i];
				index = i;
			}
		}
		if ( oczekiwany[index] == 1 )
			return true;
		else
			return false;
	}
	else
	{
		for ( int i = 0; i < dlugosc; i++ )
		{
			if ( wynik[i] > prog && oczekiwany[i] > prog )
			{

			}
			else if ( wynik[i] <= prog && oczekiwany[i] <= prog )
			{

			}
			else
			{
				return false;
			}
		}
		return true;
	}

}

bool Siec::save_to_csv ( std::string adres )
{
	std::ofstream MyExcelFile;
	MyExcelFile.open ( adres );

	if ( MyExcelFile.is_open () )
	{
		MyExcelFile << "n: " << n << " m: " << m << " l: " << l << " w: " << war_uczenia << ";\n";
		MyExcelFile << "Maciez wag warstwy wyjsciowej:;\n";
		for ( int j = 0; j < l; j++ )
		{
			MyExcelFile << j << ';';
		}
		MyExcelFile << '\n';
		for ( int i = 0; i < m; i++ )
		{
			for ( int j = 0; j < l; j++ )
			{
				MyExcelFile << wagi_war_wyj[i][j] << ";";
			}
			MyExcelFile << '\n';
		}
		MyExcelFile << "\n";
		MyExcelFile << "Biasy warstwy wyjsciowej:;\n";
		for ( int j = 0; j < m; j++ )
		{
			MyExcelFile << bias_war_wyj[j] << ';';
		}
		MyExcelFile << "\n";
		MyExcelFile << "\n";
		MyExcelFile << "Biasy warstwy Ukrytej:;\n";
		for ( int j = 0; j < l; j++ )
		{
			MyExcelFile << bias_war_ukr[j] << ';';
		}
		MyExcelFile << "\n";
		MyExcelFile << "\n";
		for ( int i = 0; i < l; i++ )
		{
			MyExcelFile << "Maciez wag warstwy ukrytej nr: " << i << "\n";
			for ( int j = 0; j < n; j++ )
			{
				MyExcelFile << wagi_war_ukr[i][j] << ";";
				if ( ( j + 1 ) % 28 == 0 ) MyExcelFile << '\n';
			}
			MyExcelFile << '\n';
		}
		return true;
	}
	else return false;

}

std::pair<double, double> Siec::badanie ( int ilosc_powotrzen, double wartosc_uczenia, double zakres_wag, double momentum, int jaka_funkcja_aktywacji, int liczba_neuronow_w_war_ukrytej )
{
	double srednia = 0;
	double wariancja = 0;
	double * his_epok = new double[ilosc_powotrzen] ();

	zmien_ilosc_neuronow_war_ukrytej ( liczba_neuronow_w_war_ukrytej );
	war_uczenia = wartosc_uczenia;
	this->momentum = momentum;

	std::streambuf* cout_sbuf = std::cout.rdbuf ();
	//std::cout.rdbuf ( NULL );
	for ( int i = 0; i < ilosc_powotrzen; i++ )
	{
		//std::cout.rdbuf ( cout_sbuf );
		std::cout << "badanie: " << i << "\n";
		//std::cout.rdbuf ( NULL );
		ustaw_wagi ( -zakres_wag, zakres_wag );
		his_epok[i] = epoki_uczenie ( 100, 3.0 );
		srednia += his_epok[i] / ilosc_powotrzen;
	}

	for ( int i = 0; i < ilosc_powotrzen; i++ )
	{
		wariancja += ( his_epok[i] - srednia )*( his_epok[i] - srednia ) / ilosc_powotrzen;
	}
	std::cout.rdbuf ( cout_sbuf );

	std::cout << srednia << " " << wariancja << "\n";

	delete[] his_epok;
	return std::make_pair ( srednia, wariancja );
}

void Siec::badanie2 ( std::string adres, int max_epok, double wartosc_uczenia, double zakres_wag, double momentum, int jaka_funkcja_aktywacji, int liczba_neuronow_w_war_ukrytej )
{
	zmien_ilosc_neuronow_war_ukrytej ( liczba_neuronow_w_war_ukrytej );
	war_uczenia = wartosc_uczenia;
	this->momentum = momentum;
	ustaw_wagi ( -zakres_wag, zakres_wag );

	std::ofstream MyExcelFile;
	std::locale mylocale ( "" );
	MyExcelFile.open ( adres );
	MyExcelFile.imbue ( mylocale );
	MyExcelFile << "epoka;suma bledow;blad walidacji;blad uczenia;blad srednio kwad walidacji;\n";

	std::cout << std::setprecision ( 4 );
	std::cout << std::fixed;
	std::cout.width ( 6 );
	double suma_bledow = 0;
	for ( int ktora_epoka = 0; ktora_epoka < max_epok; ktora_epoka++ )
	{
		suma_bledow = epoka ();

		double blad_walid = test ( ciagi_walidacyjne, ciagi_walidacyjne_war_oczekiwana, ilosc_ciagow_walidacyjnych );
		double blad_uczenia = test ( ciagi_uczace, ciagi_uczace_war_oczekiwana, ilosc_ciagow_uczacych );
		double blad_sr_walidacji = test2 ( ciagi_walidacyjne, ciagi_walidacyjne_war_oczekiwana, ilosc_ciagow_walidacyjnych );

		std::cout << "epoka: " << ktora_epoka << "\tsuma bledow prop: " << suma_bledow
			<< "\tbledne wal: " << blad_walid << "\tblad sr wal: " << blad_sr_walidacji << "\tbledne nauczone: " << blad_uczenia << "\n";

		MyExcelFile << ktora_epoka << ";" << suma_bledow << ";" << blad_walid << ";" << blad_uczenia << ";" << blad_sr_walidacji << ";\n";
	}
	MyExcelFile.close ();

}

std::pair<double*****, double*****>Siec::turbo_badanie ( int ilosc_powotrzen, std::vector<double> wartosc_uczenia, std::vector<double> zakres_wag, std::vector<double> momentum, std::vector<int> jaka_funkcja_aktywacji, std::vector<int> liczba_neuronow_w_war_ukrytej )
{
	std::streambuf* cout_sbuf = std::cout.rdbuf ();
	std::cout << "max: " << liczba_neuronow_w_war_ukrytej.size ()* jaka_funkcja_aktywacji.size () * momentum.size ()* zakres_wag.size ()*wartosc_uczenia.size () << "\n";

	double***** wyniki_epoki = new double ****[(int) wartosc_uczenia.size ()];
	double***** wyniki_wariancje = new double ****[(int) wartosc_uczenia.size ()];
	std::pair <double, double> wynik;

	for ( int ii = 0; ii < (int) wartosc_uczenia.size (); ii++ )
	{
		wyniki_epoki[ii] = new double ***[(int) zakres_wag.size ()];
		wyniki_wariancje[ii] = new double ***[(int) zakres_wag.size ()];
		for ( int jj = 0; jj < (int) zakres_wag.size (); jj++ )
		{
			wyniki_epoki[ii][jj] = new double **[(int) momentum.size ()];
			wyniki_wariancje[ii][jj] = new double **[(int) momentum.size ()];
			for ( int kk = 0; kk < (int) momentum.size (); kk++ )
			{
				wyniki_epoki[ii][jj][kk] = new double *[(int) jaka_funkcja_aktywacji.size ()];
				wyniki_wariancje[ii][jj][kk] = new double *[(int) jaka_funkcja_aktywacji.size ()];
				for ( int ll = 0; ll < (int) jaka_funkcja_aktywacji.size (); ll++ )
				{
					wyniki_epoki[ii][jj][kk][ll] = new double[(int) liczba_neuronow_w_war_ukrytej.size ()];
					wyniki_wariancje[ii][jj][kk][ll] = new double[(int) liczba_neuronow_w_war_ukrytej.size ()];
					for ( int mm = 0; mm < (int) liczba_neuronow_w_war_ukrytej.size (); mm++ )
					{
						std::cout << (int) mm +
							liczba_neuronow_w_war_ukrytej.size ()*ll +
							liczba_neuronow_w_war_ukrytej.size ()* jaka_funkcja_aktywacji.size () * kk +
							liczba_neuronow_w_war_ukrytej.size ()* jaka_funkcja_aktywacji.size () * momentum.size ()* jj +
							liczba_neuronow_w_war_ukrytej.size ()* jaka_funkcja_aktywacji.size () * momentum.size ()* zakres_wag.size ()*ii << " ";

						std::cout.rdbuf ( NULL );
						wynik = badanie ( ilosc_powotrzen, wartosc_uczenia[ii], zakres_wag[jj], momentum[kk], jaka_funkcja_aktywacji[ll], liczba_neuronow_w_war_ukrytej[mm] );
						wyniki_epoki[ii][jj][kk][ll][mm] = wynik.first;
						wyniki_wariancje[ii][jj][kk][ll][mm] = wynik.second;
						std::cout.rdbuf ( cout_sbuf );
					}
				}
			}
		}
	}
	std::cout.rdbuf ( cout_sbuf );
	return std::make_pair ( wyniki_epoki, wyniki_wariancje );
}


void Siec::dawaj_skal_wagi_ukr ( std::string adres )
{
	//operator pliku
	std::ofstream MyExcelFile;
	MyExcelFile.open ( adres );

	//przechodzimy po wszystkich neuronach w warstwie ukrytej
	for ( int i = 0; i < l; i++ )
	{
		double max = *std::min_element ( wagi_war_ukr[i], wagi_war_ukr[i] + m );// znajdujemy max wage wychodzaca z tego neuronu
		double min = *std::max_element ( wagi_war_ukr[i], wagi_war_ukr[i] + m );// znajdujemy min wage wychodzaca z tego neuronu

		// w excelu komorka po komoce robimy z wektora wag[784] maciez[28][28]
		for ( int j = 0; j < 28; j++ )
		{
			for ( int k = 0; k < 28; k++ )
			{
				/*                                                   \/ skalujemy w dol aby najwieksza waga byla 1, a najmniejsza 0 */
				MyExcelFile << ( wagi_war_ukr[i][j * 28 + k] - min ) / ( max - min ) << ";";
			}
			MyExcelFile << "\n";
		}
		MyExcelFile << "\n";
	}

	MyExcelFile.close ();
}

void Siec::dawaj_skal_wag_wyj ( std::string adres )
{
	std::ofstream MyExcelFile;
	MyExcelFile.open ( adres );
	for ( int i = 0; i < l; i++ )
	{
		double max = 1;
		double min = 0;
		for ( int j = 0; j < 28; j++ )
		{
			for ( int k = 0; k < 28; k++ )
			{
				MyExcelFile << ( wagi_war_wyj[j * 28 + k][i] - min ) / ( max - min ) << ";";
			}
			MyExcelFile << "\n";
		}
		MyExcelFile << "\n";
	}

	MyExcelFile.close ();
}


/*
notatki z msid
Ciag treningowy - innaczej ci¹g ucz¹cy, na którego podstawie liczone s¹ parametry modelu
Ci¹g walidacyjny - s³u¿y do porównania modeli (róznych stopni wielomianów)
Ci¹g testowy - s³u¿y do sprawdzenia dzia³ania modelu
*/

std::ostream & operator<<( std::ostream & Strm_Wyj, const Siec & siec )
{
	Strm_Wyj << "n: " << siec.n << " m: " << siec.m << " l: " << siec.l << " w: " << siec.war_uczenia << "\n";
	Strm_Wyj << "Maciez wag warstwy wyjsciowej:\n";
	for ( int i = 0; i < siec.m; i++ )
	{
		for ( int j = 0; j < siec.l; j++ )
		{
			Strm_Wyj << siec.wagi_war_wyj[i][j] << " ";
		}
		Strm_Wyj << '\n';
	}
	Strm_Wyj << "Maciez wag warstwy ukrytej:\n";
	for ( int i = 0; i < siec.l; i++ )
	{
		for ( int j = 0; j < siec.n; j++ )
		{
			Strm_Wyj << siec.wagi_war_ukr[i][j] << " ";
		}
		Strm_Wyj << '\n';
	}
	Strm_Wyj << "Maciez biasow warstwy wyjsciowej:\n";
	for ( int i = 0; i < siec.m; i++ )
	{
		Strm_Wyj << siec.bias_war_wyj[i] << " ";
	}
	Strm_Wyj << '\n';
	Strm_Wyj << "Maciez biasow warstwy ukrytej:\n";
	for ( int i = 0; i < siec.l; i++ )
	{
		Strm_Wyj << siec.bias_war_ukr[i] << " ";
	}
	return Strm_Wyj;
}

std::istream & operator >> ( std::istream & Strm_We, Siec & siec )
{
	int n, m, l;
	double war_ucz;
	Strm_We.ignore ( 3 );
	Strm_We >> n;
	Strm_We.ignore ( 3 );
	Strm_We >> m;
	Strm_We.ignore ( 3 );
	Strm_We >> l;
	Strm_We.ignore ( 3 );
	Strm_We >> war_ucz;
	Strm_We.ignore ( 32 );

	siec = Siec ( n, m, l, war_ucz );

	for ( int i = 0; i < m && !Strm_We.eof (); i++ )
	{
		for ( int j = 0; j < l && !Strm_We.eof (); j++ )
		{
			Strm_We >> siec.wagi_war_wyj[i][j];
		}
		Strm_We.ignore ( 1 );
	}

	Strm_We.ignore ( 28 );
	for ( int i = 0; i < l && !Strm_We.eof (); i++ )
	{
		for ( int j = 0; j < n && !Strm_We.eof (); j++ )
		{
			Strm_We >> siec.wagi_war_ukr[i][j];
		}
		Strm_We.ignore ( 1 );
	}
	Strm_We.ignore ( 34 );
	for ( int i = 0; i < m && !Strm_We.eof (); i++ )
	{
		Strm_We >> siec.bias_war_wyj[i];
	}
	Strm_We.ignore ( 1 );

	Strm_We.ignore ( 31 );
	for ( int i = 0; i < l && !Strm_We.eof (); i++ )
	{
		Strm_We >> siec.bias_war_ukr[i];
	}

	return Strm_We;
}
