/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Colloid transport

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>

struct UserData
{

    double          Kn; 
    double          Gn;
    Vec3_t           g;
    Array<double>   Fb;
    std::ofstream      oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0;i<dom.Disks.Size();i++)
    {
        dom.Disks[i]->Ff = dom.Disks[i]->M*dat.g;
        double delta;
        //delta =   dat.Xmin(0) - dom.Disks[i]->X(0) + dom.Disks[i]->R;
        //if (delta > 0.0)  dom.Disks[i]->Ff(0) += dat.Kn*delta;
        delta = dom.Disks[i]->R - dom.Disks[i]->X(1);
        if (delta > 0.0)
        {
            dom.Disks[i]->Ff(1) += dat.Kn*delta-dat.Gn*dom.Disks[i]->V(1);
            dat.Fb[i] += dom.Disks[i]->V(0)*dom.dt;
            if (fabs(dat.Fb[i])>=dom.Disks[i]->Mu*delta) dat.Fb[i] *= dom.Disks[i]->Mu*delta/fabs(dat.Fb[i]);
            dom.Disks[i]->Ff(0) += -dat.Kn*dat.Fb[i];
        }
    }
}

void Report(LBM::Domain & dom, void * UD)
{
    //std::cout << dom.Disks[0]->X << dom.Disks[1]->X << std::endl;        
}

int main(int argc, char **argv) try
{
    size_t Nproc  = 1; 
    size_t nx     = 200;
    size_t ny     = 1000;
    double nu     = 1.5e-5;
    double dx     = 4.0e-6;
    double dt     = 2.0e-7;    
    double rc     = 2.5e-5;
    double rho    = 2500.0;
    double rhof   = 1000.0;
    double Kn     = 0.01*M_PI*rho*rc*rc/(dt*dt);
    double Tf     = 0.05;
    if (argc>=2) Nproc = atoi(argv[1]);
    if (argc>=3) Tf    = atof(argv[2]);

    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Alpha    = 0.1*rc;
    //Dom.RotPar   = false;
    //Dom.Sc       = 0.0;
    dat.g        = 0.0,-9.8,0.0;
    dat.Kn       = Kn;
    dat.Gn       = -0.2;
    double me    = M_PI*rho*rc*rc;
    dat.Gn = 2.0*sqrt((pow(log(-dat.Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-dat.Gn),2.0)))*me;


    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }

    Dom.AddDisk(-1,Vec3_t(0.5*nx*dx-2.0*rc,0.2*ny*dx,0.0),OrthoSys::O,OrthoSys::O,rho,rc,dt);
    Dom.AddDisk(-1,Vec3_t(0.5*nx*dx+2.0*rc,0.2*ny*dx,0.0),OrthoSys::O,OrthoSys::O,rho,rc,dt);
    //Dom.AddDisk(-1,Vec3_t(0.5*nx*dx-2.0*rc,rc,0.0),OrthoSys::O,OrthoSys::O,rho,rc,dt);
    //Dom.AddDisk(-1,Vec3_t(0.5*nx*dx+2.0*rc,rc,0.0),OrthoSys::O,OrthoSys::O,rho,rc,dt);

    for (size_t i=0;i<Dom.Disks.Size();i++)
    {
        Dom.Disks[i]->Kn = Kn;
        Dom.Disks[i]->Kt = Kn;
        Dom.Disks[i]->Gn = -0.5;
        Dom.Disks[i]->Mu =  0.5;
    }

    dat.Fb.Resize(Dom.Disks.Size());
    

    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rhof, OrthoSys::O);
    }

    //Solving
    Dom.Solve(Tf,0.01*Tf,Setup,Report,"column",true,Nproc);
    dat.oss_ss.close();
 
}
MECHSYS_CATCH

