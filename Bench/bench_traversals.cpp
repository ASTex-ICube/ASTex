/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



#include <iostream>
#include <chrono>

#include <ASTex/image_gray.h>

using namespace ASTex;


const int NB_ITER = 40;
const int IMG_SZ = 2048;


struct PipoImg
{
	int width;
	int height;
	std::vector<double> data;

	PipoImg(int w, int h):
		width(w),
		height(h),
		data(w*h)
	{
	}

	inline double& pix(int i, int j) { return data[j*width+i];}
};

void bench_table()
{
	PipoImg im(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	for(int j=0; j<im.height; ++j)
		for(int i=0; i<im.width; ++i)
			im.pix(i,j)=0.0;

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0; j<im.height; ++j)
			for(int i=0; i<im.width; ++i)
				im.pix(i,j)+=0.1;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0; j<im.height; ++j)
			for(int i=0; i<im.width; ++i)
				im.pix(i,j)-=0.1;
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "table[i,j] : "<< elapsed_seconds.count() << " s."<< std::endl;

	start_chrono = std::chrono::system_clock::now();

	for(double& P: im.data)
		P=0.0;

	for (int k=0;k<NB_ITER;++k)
	{
		for(double& P: im.data)
			P=+0.1;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for(double& P: im.data)
			P=-0.1;
	}

	elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "table[i] : "<< elapsed_seconds.count() << " s."<< std::endl;

	start_chrono = std::chrono::system_clock::now();

	std::vector<double>::iterator P = im.data.begin();
	std::vector<double>::iterator Q = im.data.end();

	while (P!=Q)
	{
		*P++ = 0;
	}


	for (int k=0;k<NB_ITER;++k)
	{
		P = im.data.begin();
		while (P!=Q)
		{
			*P++ += 0.1;
		}
	}

	for (int k=0;k<NB_ITER;++k)
	{
		P = im.data.begin();
		while (P!=Q)
		{
			*P++ -= 0.1;
		}
	}

	elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "*p++ : "<< elapsed_seconds.count() << " s."<< std::endl;

}




void bench_pixel_ptr()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	double* Pi = img.getPixelsPtr();
	double* Q = Pi+(img.width()*img.height());

	double* P = Pi;
	while (P!=Q)
	{
		*P++ = 0;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		P = Pi;
		while (P!=Q)
		{
			*P++ += 0.1;
		}
	}

	for (int k=0;k<NB_ITER;++k)
	{
		P = Pi;
		while (P!=Q)
		{
			*P++ -= 0.1;
		}
	}


	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "pixel_ptr : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});
}


void bench_pixel_absolute()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	for(int j=0; j<img.height(); ++j)
		for(int i=0; i<img.width(); ++i)
			img.pixelAbsolute(i,j)=0.0;

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0; j<img.height(); ++j)
			for(int i=0; i<img.width(); ++i)
				img.pixelAbsolute(i,j) += 0.1;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0; j<img.height(); ++j)
			for(int i=0; i<img.width(); ++i)
				img.pixelAbsolute(i,j) -= 0.1;
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "pixel_absolute : "<< elapsed_seconds.count() << " s."<< std::endl;


	start_chrono = std::chrono::system_clock::now();

	for(int j=0,h=img.height(); j<h; ++j)
		for(int i=0,w=img.width(); i<w; ++i)
			img.pixelAbsolute(i,j)=0.0;

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0,h=img.height(); j<h; ++j)
			for(int i=0,w=img.width(); i<w; ++i)
				img.pixelAbsolute(i,j)+=0.1;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0,h=img.height(); j<h; ++j)
			for(int i=0,w=img.width(); i<w; ++i)
				img.pixelAbsolute(i,j)-=0.1;
	}

	elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "pixel_absolute opt: "<< elapsed_seconds.count() << " s."<< std::endl;



	start_chrono = std::chrono::system_clock::now();

	for_indices(0,img.width(),0,img.height(),[&](int i, int j)
	{
		img.pixelAbsolute(i,j)=0.0;
	});

	for (int k=0;k<NB_ITER;++k)
	{
		for_indices(0,img.width(),0,img.height(),[&](int i, int j)
		{
			img.pixelAbsolute(i,j) += 0.1;
		});
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for_indices(0,img.width(),0,img.height(),[&](int i, int j)
		{
			img.pixelAbsolute(i,j) -= 0.1;
		});
	}

	elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "pixel_absolute indices : "<< elapsed_seconds.count() << " s."<< std::endl;




	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});
}


void bench_pixel_absolute_xy()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	for(int i=0,w=img.width(); i<w; ++i)
		for(int j=0,h=img.height(); j<h; ++j)
			img.pixelAbsolute(i,j)=0.0;

	for (int k=0;k<NB_ITER;++k)
	{
		for(int i=0,w=img.width(); i<w; ++i)
			for(int j=0,h=img.height(); j<h; ++j)
				img.pixelAbsolute(i,j)+=0.1;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for(int i=0,w=img.width(); i<w; ++i)
			for(int j=0,h=img.height(); j<h; ++j)
				img.pixelAbsolute(i,j)-=0.1;
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "pixel_absolute xy : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});
}



void bench_pixel_relative()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	for(int j=0,h=img.height(); j<h; ++j)
		for(int i=0,w=img.width(); i<w; ++i)
			img.pixelRelative(i,j)=0.0;

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0,h=img.height(); j<h; ++j)
			for(int i=0,w=img.width(); i<w; ++i)
				img.pixelRelative(i,j)+=0.1;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for(int j=0,h=img.height(); j<h; ++j)
			for(int i=0,w=img.width(); i<w; ++i)
				img.pixelRelative(i,j)-=0.1;
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "pixel_relative : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});
}


void bench_iterator()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	for (auto it = img.beginIterator(); !it.IsAtEnd(); ++it)
	{
		it.Value() = 0.0;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for (auto it = img.beginIterator(); !it.IsAtEnd(); ++it)
		{
			it.Value() += 0.1;
		}
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for (auto it = img.beginIterator(); !it.IsAtEnd(); ++it)
		{
			it.Value() -= 0.1;
		}
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "fast_iterator : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});
}


void bench_rev_iterator()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	for (auto it = img.endIterator(); !it.IsAtBegin();)
	{
		--it;
		it.Value() = 0.0;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for (auto it = img.endIterator(); !it.IsAtBegin();)
		{
			--it;
			it.Value() += 0.1;
		}
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for (auto it = img.endIterator(); !it.IsAtBegin();)
		{
			--it;
			it.Value() -= 0.1;
		}
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "fast_iterator backward : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});
}

void bench_iterator_indexed()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	for (auto it = img.beginIteratorIndexed(); !it.IsAtEnd(); ++it)
	{
		it.Value() = 0.0;
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for (auto it = img.beginIteratorIndexed(); !it.IsAtEnd(); ++it)
		{
			it.Value() += 0.1;
		}
	}

	for (int k=0;k<NB_ITER;++k)
	{
		for (auto it = img.beginIteratorIndexed(); !it.IsAtEnd(); ++it)
		{
			it.Value() -= 0.1;
		}
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "iterator indexed : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});

}

void bench_for_all_pixels()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);


	auto start_chrono = std::chrono::system_clock::now();

	img.for_all_pixels([] (double& p)
	{
		p = 0.0;
	});


	for (int k=0;k<NB_ITER;++k)
	{
		img.for_all_pixels([] (double& p)
		{
			p += 0.1;
		});
	}

	for (int k=0;k<NB_ITER;++k)
	{
		img.for_all_pixels([] (double& p)
		{
			p -= 0.1;
		});
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "for_all_pixels : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});

}


void bench_for_all_pixels_mask()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);


	auto start_chrono = std::chrono::system_clock::now();

	img.for_all_pixels([] (double& p)
	{
		p = 0.0;
	},[] (int x,int y) { return (x>=0) && (y>=0);});


	for (int k=0;k<NB_ITER;++k)
	{
		img.for_all_pixels([] (double& p)
		{
			p += 0.1;
		},[] (int x,int y) { return (x>=0) && (y>=0);});
	}

	for (int k=0;k<NB_ITER;++k)
	{
		img.for_all_pixels([] (double& p)
		{
			p -= 0.1;
		},[] (int x,int y) { return (x>=0) && (y>=0);});
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "for_all_pixels_mask : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
			std::cout << "pb p = "<<p<< std::endl;
	});

}



void bench_parallel_for_all_pixels()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	/*std::chrono::time_point<std::chrono::system_clock>*/
	auto start_chrono = std::chrono::system_clock::now();

	img.parallel_for_all_pixels([] (double& p)
	{
		p = 0.0;
	});


	for (int k=0;k<NB_ITER;++k)
	{
		img.parallel_for_all_pixels([] (double& p)
		{
			p += 0.1;
		});
	}

	for (int k=0;k<NB_ITER;++k)
	{
		img.parallel_for_all_pixels([] (double& p)
		{
			p -= 0.1;
		});
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "parallel_for_all_pixels : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
		{
			std::cout << "pb parallel_for_all_pixels p = "<<p<< std::endl;
		}
	});

}



void bench_parallel_for_all_pixels_mask()
{
	ImageGrayd img(IMG_SZ,IMG_SZ);

	auto start_chrono = std::chrono::system_clock::now();

	img.parallel_for_all_pixels([] (double& p)
	{
		p = 0.0;
	},[] (int x,int y) { return (x>=0) && (y>=0);});


	for (int k=0;k<NB_ITER;++k)
	{
		img.parallel_for_all_pixels([] (double& p)
		{
			p += 0.1;
		},[] (int x,int y) { return (x>=0) && (y>=0);});
	}

	for (int k=0;k<NB_ITER;++k)
	{
		img.parallel_for_all_pixels([] (double& p)
		{
			p -= 0.1;
		},[] (int x,int y) { return (x>=0) && (y>=0);});
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "parallel_for_all_pixels_mask : "<< elapsed_seconds.count() << " s."<< std::endl;

	img.for_all_pixels([] (double p)
	{
		if (std::abs(p) > 0.0000000001)
		{
			std::cout << "pb parallel_for_all_pixels p = "<<p<< std::endl;
		}
	});

}


int main()
{
	bench_table();

	bench_pixel_ptr();

	bench_pixel_absolute();

	bench_pixel_absolute_xy();

	bench_pixel_relative();

	bench_iterator();

	bench_rev_iterator();

	bench_iterator_indexed();

	bench_for_all_pixels();

	bench_for_all_pixels_mask();

	bench_parallel_for_all_pixels();

	bench_parallel_for_all_pixels_mask();

	return EXIT_SUCCESS;
}

