#include "Tools.h"


/*int z zakresu <min, max)*/
int Tools::r_int ( int min, int max )
{
	int f = rand () % max + min;
	return f;
}

/*int z zakresu <0, 32767>*/
int Tools::r_int ()
{
	return rand ();
}

/*double z zakresu <min,max>*/
double Tools::r_db ( double min, double max )
{
	double f = (double) rand () / RAND_MAX;
	return min + f * ( max - min );
}

/*double z zakresu <0,1> */
double Tools::r_db ()
{
	return (double) rand () / RAND_MAX;
}

/*Inicjalizacja losowych liczb*/
bool Tools::init ()
{
	srand ( (unsigned int) time ( NULL ) );
	return true;
}



std::vector<std::string> Tools::get_all_files_names_within_folder ( std::string folder )
{
	std::vector<std::string> names;
	std::string search_path = folder + "/*.*";
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile ( search_path.c_str (), &fd );
	if ( hFind != INVALID_HANDLE_VALUE )
	{
		do
		{
			// read all (real) files in current folder
			// , delete '!' read other 2 default folder . and ..
			if ( !( fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY ) )
			{
				names.push_back ( fd.cFileName );
			}
		} while ( ::FindNextFile ( hFind, &fd ) );
		::FindClose ( hFind );
	}
	return names;
}

/*srednie rozwiazanie, duza zlozonosc obliczeniowa*/
std::vector<int> Tools::vektor_roznych_liczb ( int min, int max, int ilosc )
{
	if ( max - min < ilosc )
	{
		std::cout << " Tools::vektor_roznych_liczb - error: wymagana ilosc wieksza niz mozliwy zakres\n";
		throw std::exception ( "wymagana ilosc wieksza niz mozliwy zakres" );
	}

	static std::vector<int> generatedValues;

	for ( int i = 0; i < ilosc; i++ )
	{
		std::cout << "for in tools vektor\n";
		int num = rand () % max + min;

		bool contains = false;
		for ( auto & p : generatedValues )
		{
			if ( p == num ) contains = true;
		}
		while ( contains )
		{
			std::cout << "while  in tools vektor\n";
			num = rand () % max + min;
			contains = false;
			for ( auto & p : generatedValues )
			{
				if ( p == num ) contains = true;
			}
		}
		generatedValues.push_back ( num );
	}
	return generatedValues;
}

//zlozonosc obliczeniowa O( 2(max-min) )
std::vector<int> Tools::vektor_roznych_liczb_ver2 ( int min, int max )
{
	std::vector<int> generatedValues;
	for ( int i = min; i < max; i++ )
	{
		generatedValues.push_back ( i );
	}
	std::random_shuffle ( generatedValues.begin (), generatedValues.end () );
	return generatedValues;
}

int * Tools::img_to_array ( std::string adres_katalogu, std::string nazwa_pliku, int & dlugosc )
{
	//std::cout << adres_katalogu + nazwa_pliku << '\n';
	cv::Mat image = cv::imread ( adres_katalogu + nazwa_pliku );
	dlugosc = ( image.rows*image.cols ) + 1;
	//std::cout << "" << image.rows << " " << image.cols <<" " << ( image.rows*image.cols ) + 1<<  '\n';
	int * tab = new int[dlugosc] ();
	tab[0] = nazwa_pliku[0] - 48;
	for ( int i = 0; i < image.rows; i++ )
	{
		for ( int j = 0; j < image.cols; j++ )
		{
			cv::Vec3b bgrPixel = image.at<cv::Vec3b> ( i, j );
			if ( bgrPixel[0] < 100 )
				tab[1 + image.cols*i + j] = 1;
			else
				tab[1 + image.cols*i + j] = 0;
		}
	}
	return tab;
}

int ** Tools::folder_of_img_to_array ( std::string adres_katalogu, int & ilosc_w_folderze, int & dlugosc_vektorow )
{
	std::vector <std::string> vektor_nazw = Tools::get_all_files_names_within_folder ( adres_katalogu );
	std::cout << adres_katalogu << " " << vektor_nazw.size () << '\n';
	ilosc_w_folderze = (int) vektor_nazw.size ();

	//przelosowanie ciagow
	//std::random_shuffle ( vektor_nazw.begin (), vektor_nazw.end () );

	int ** tab = new int *[ilosc_w_folderze];
	int i = 0;
	for ( auto & p : vektor_nazw )
	{
		//tab[i] = new int[71]; // powino byc innaczej niz podane stale ale chwilowo najlepsze rozwiazanie 
		tab[i] = img_to_array ( adres_katalogu + "/", p, dlugosc_vektorow );
		i++;
	}
	std::cout << "Wczytywanie z folderu zakonczone\n";
	return tab;
}

std::string Tools::int_to_fancy_vector ( int i )
{
	switch ( i )
	{
	case 0:
		return std::string ( "1;0;0;0;0;0;0;0;0;0;" );
	case 1:
		return std::string ( "0;1;0;0;0;0;0;0;0;0;" );
	case 2:
		return std::string ( "0;0;1;0;0;0;0;0;0;0;" );
	case 3:
		return std::string ( "0;0;0;1;0;0;0;0;0;0;" );
	case 4:
		return std::string ( "0;0;0;0;1;0;0;0;0;0;" );
	case 5:
		return std::string ( "0;0;0;0;0;1;0;0;0;0;" );
	case 6:
		return std::string ( "0;0;0;0;0;0;1;0;0;0;" );
	case 7:
		return std::string ( "0;0;0;0;0;0;0;1;0;0;" );
	case 8:
		return std::string ( "0;0;0;0;0;0;0;0;1;0;" );
	case 9:
		return std::string ( "0;0;0;0;0;0;0;0;0;1;" );
	default:
		return std::string ();
	}
}

int Tools::reverseInt ( int i )
{
	unsigned char c1, c2, c3, c4;

	c1 = i & 255;
	c2 = ( i >> 8 ) & 255;
	c3 = ( i >> 16 ) & 255;
	c4 = ( i >> 24 ) & 255;

	return ( (int) c1 << 24 ) + ( (int) c2 << 16 ) + ( (int) c3 << 8 ) + c4;
}

void Tools::read_mnist_data ( std::string adres, int & ilosc_obrazow, int & wielkosc_obrazow, double ** & data )
{
	std::ifstream file ( adres, std::ios::binary );
	if ( file.is_open () )
	{
		int magic_number = 0;
		ilosc_obrazow = 0;
		int n_rows = 0;
		int n_cols = 0;

		file.read ( (char*) &magic_number, sizeof ( magic_number ) );
		magic_number = reverseInt ( magic_number );
		file.read ( (char*) &ilosc_obrazow, sizeof ( ilosc_obrazow ) );
		ilosc_obrazow = reverseInt ( ilosc_obrazow );
		file.read ( (char*) &n_rows, sizeof ( n_rows ) );
		n_rows = reverseInt ( n_rows );
		file.read ( (char*) &n_cols, sizeof ( n_cols ) );
		n_cols = reverseInt ( n_cols );

		wielkosc_obrazow = n_rows * n_cols;

		data = new double *[ilosc_obrazow];
		for ( int i = 0; i < ilosc_obrazow; i++ )
		{
			data[i] = new double[wielkosc_obrazow];
			for ( int r = 0; r < n_rows; ++r )
			{
				for ( int c = 0; c < n_cols; ++c )
				{
					unsigned char temp = 0;
					file.read ( (char*) &temp, 1 );
					data[i][( n_cols*r ) + c] = (double) temp / 255;
				}
			}
		}
	}
	else
	{
		std::cout << "read_mnist: File is not open\n";
	}
}

void Tools::read_mnist_labels ( std::string adres, int & ilosc_etykiet, double **& etykiety )
{

	std::ifstream file ( adres, std::ios::binary );
	if ( file.is_open () )
	{
		int magic_number = 0;
		ilosc_etykiet = 0;

		file.read ( (char*) &magic_number, sizeof ( magic_number ) );
		magic_number = reverseInt ( magic_number );
		file.read ( (char*) & ilosc_etykiet, sizeof ( ilosc_etykiet ) );
		ilosc_etykiet = reverseInt ( ilosc_etykiet );

		etykiety = new double*[ilosc_etykiet];
		for ( int i = 0; i < ilosc_etykiet; i++ )
		{
			etykiety[i] = new double[10]();
			unsigned char temp = 0;
			file.read ( (char*) &temp, 1 );
			etykiety[i][(int) temp] = 1;
		}
	}
	else
	{
		std::cout << "read_mnist_label: File is not open\n";
	}
}
void Tools::wagi_to_png ( std::string nazwa, int ile_wag, int wys, int szer, double ** wagi )
{
	int img_size = (int) ceil ( std::sqrt ( ile_wag ) );
	int img_hi = img_size*wys + ( img_size - 1 ) * 4;
	int img_wid = img_size*szer + ( img_size - 1 ) * 4;

	cv::Mat_<cv::Vec3b> img ( img_hi, img_wid, cv::Vec3b ( 125, 0, 0 ) );
	int iter_wid = 0;
	int iter_hi = 0;

	for ( int i = 0; i < ile_wag; i++ )
	{
		double max = *std::min_element ( wagi[i], wagi[i] + wys*szer );
		double min = *std::max_element ( wagi[i], wagi[i] + wys*szer );
		for ( int j = 0; j < wys; j++ )
		{
			for ( int k = 0; k < szer; k++ )
			{
				int kolor = (int) (255 - ( wagi[i][j * 28 + k] - min ) / ( max - min ) * 255);
				img ( iter_hi*( wys + 4 ) + j, iter_wid*( szer + 4 ) + k )[0] = kolor;
				img ( iter_hi*( wys + 4 ) + j, iter_wid*( szer + 4 ) + k )[1] = kolor;
				img ( iter_hi*( wys + 4 ) + j, iter_wid*( szer + 4 ) + k )[2] = kolor;
			}
		}
		iter_wid++;
		if ( iter_wid == img_size )
		{
			iter_wid = 0;
			iter_hi++;
		}
	}
	try
	{
		cv::imwrite ( nazwa, img );
	}
	catch ( std::runtime_error& ex )
	{
		fprintf ( stderr, "Exception converting image to PNG format: %s\n", ex.what () );
	}
}

void Tools::wagi_to_png_transp ( std::string nazwa, int ile_wag, int wys, int szer, double ** wagi )
{
	double ** tab = new double *[ile_wag];
	for ( int i = 0; i < ile_wag; i++ )
	{
		tab[i] = new double[wys*szer];
		for ( int j = 0; j < wys*szer; j++ )
		{
			tab[i][j] = wagi[j][i];
		}
	}
	wagi_to_png ( nazwa, ile_wag, wys, szer, tab );
	delete[] tab;
}


/*
	zlozonosc obliczeniowa O( 3(max-min) ) ale zwraca tylko vektor o podanej dlugosci z zakresu
*/
std::vector<int> Tools::vektor_roznych_liczb_ver2 ( int min, int max, int ilosc )
{
	std::vector<int> generatedValues = vektor_roznych_liczb_ver2 ( min, max );
	generatedValues.resize ( ilosc );
	return generatedValues;
}

bool Tools::static_init ()
{
	srand ( (unsigned int) time ( NULL ) );
	return true;
}
