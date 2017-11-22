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
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace ASTex;


const int NB_ITER = 400;

const int NB=2000000;

void bench_normal()
{
	std::vector<Eigen::Vector3d> Va(NB);
	std::vector<Eigen::Vector3d> Vb(NB);
	std::vector<Eigen::Vector3d> Vc(NB);

	auto start_chrono = std::chrono::system_clock::now();

	for(int i=0; i<NB; ++i)
		Va[i]=Eigen::Vector3d(1.77,2.55,3.11);
	for(int i=0; i<NB; ++i)
		Vb[i]=Eigen::Vector3d(7.7,-5.5,1.12);
	for(int j=0; j<NB_ITER; ++j)
	{
		for(int i=0; i<NB; ++i)
			Vc[i] = 3.0*Vc[i] +(2*Va[i].cross(4*Vb[i]));
		Vc.swap(Vb);
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "normal : "<< elapsed_seconds.count() << " s."<< std::endl;
}

void bench_omp()
{
	std::vector<Eigen::Vector3d> Va(NB);
	std::vector<Eigen::Vector3d> Vb(NB);
	std::vector<Eigen::Vector3d> Vc(NB);

	auto start_chrono = std::chrono::system_clock::now();

#pragma omp parallel for
	for(int i=0; i<NB; ++i)
		Va[i]=Eigen::Vector3d(1.77,2.55,3.11);
#pragma omp parallel for
	for(int i=0; i<NB; ++i)
		Vb[i]=Eigen::Vector3d(7.7,-5.5,1.12);
	for(int j=0; j<NB_ITER; ++j)
	{
#pragma omp parallel for
		for(int i=0; i<NB; ++i)
			Vc[i] = 3.0*Vc[i] +(2*Va[i].cross(4*Vb[i]));
		Vc.swap(Vb);
	}

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "omp : "<< elapsed_seconds.count() << " s."<< std::endl;
}

void bench_parallel_for()
{
	std::vector<Eigen::Vector3d> Va(NB);
	std::vector<Eigen::Vector3d> Vb(NB);
	std::vector<Eigen::Vector3d> Vc(NB);

	auto start_chrono = std::chrono::system_clock::now();

	ASTex::parallel_for(0,NB,[&](int i)
	{
		Va[i]=Eigen::Vector3d(1.77,2.55,3.11);
	});
	ASTex::parallel_for(0,NB,[&](int i)
	{
		Vb[i]=Eigen::Vector3d(7.7,-5.5,1.12);
	});
	for(int j=0; j<NB_ITER; ++j)
	{
		ASTex::parallel_for(0,NB,[&](int i)
		{
			Vc[i] = 3.0*Vc[i] +(2*Va[i].cross(4*Vb[i]));
		});
		Vc.swap(Vb);
	}


	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "parallel_for : "<< elapsed_seconds.count() << " s."<< std::endl;
}



int main()
{
	bench_normal();

	bench_omp();

	bench_parallel_for();

	return EXIT_SUCCESS;
}

