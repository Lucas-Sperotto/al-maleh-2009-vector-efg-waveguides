// SPDX-License-Identifier: MIT
#include "al_maleh/spectral.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using namespace al_maleh;
namespace {
void write_field(const std::filesystem::path& file,
                 const MLSEvaluator& mls,
                 const Eigen::VectorXd& v, Rectangle rect,
                 int nx,int ny){
  std::ofstream stream(file);
  if(!stream)throw std::runtime_error("Cannot create field CSV");
  stream<<std::setprecision(17)<<"x,y,Hx,Hy,divH,curlH\n";
  const auto N=static_cast<Eigen::Index>(mls.nodes().size());
  for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){
    const Point2 p{rect.xmin+(i+.5)*(rect.xmax-rect.xmin)/nx,
                   rect.ymin+(j+.5)*(rect.ymax-rect.ymin)/ny};
    const auto s=mls.evaluate(p);
    double hx=0,hy=0,div=0,curl=0;
    for(std::size_t k=0;k<s.active_nodes.size();++k){
      const auto id=static_cast<Eigen::Index>(s.active_nodes[k]);
      hx+=s.phi[k]*v[id];hy+=s.phi[k]*v[N+id];
      div+=s.dx[k]*v[id]+s.dy[k]*v[N+id];
      curl+=s.dx[k]*v[N+id]-s.dy[k]*v[id];
    }
    stream<<p.x<<','<<p.y<<','<<hx<<','<<hy<<','
          <<div<<','<<curl<<'\n';
  }
}
}
int main(int argc,char** argv){
 try {
   // cartesian_spectrum [nx ny [degree [support_factor [output_dir]]]]
   if(argc>6)throw std::invalid_argument("Too many arguments");
   const int nx=argc>1?std::stoi(argv[1]):5;
   const int ny=argc>2?std::stoi(argv[2]):5;
   const int degree=argc>3?std::stoi(argv[3]):1;
   const double support=argc>4?std::stod(argv[4]):3.3;
   const std::filesystem::path out=argc>5?argv[5]:"out/m3/rect";
   if(nx<3||ny<3||nx>17||ny>17||(degree!=1&&degree!=2)||
      !std::isfinite(support)||support<=0)
     throw std::invalid_argument(
       "Usage: cartesian_spectrum [nx>=3 ny>=3 degree=1|2 support_factor>0 [out_dir]]");
   const Rectangle rect{0,2,0,1};
   const double h=std::max(2.0/(nx-1),1.0/(ny-1));
   const double radius=support*h;
   const MLSEvaluator mls(make_regular_cloud(rect,nx,ny,radius),degree);
   const BackgroundQuadrature quad{rect,nx-1,ny-1,3};
   const auto a=assemble_cartesian(mls,quad);
   const auto s=solve_generalized_symmetric(a.stiffness,a.mass);
   std::filesystem::create_directories(out);
   std::ofstream modes(out/"eigenvalues.csv");
   if(!modes)throw std::runtime_error("Cannot create eigenvalues.csv");
   modes<<std::setprecision(17)<<"index,lambda,kc,status,relative_residual\n";
   for(Eigen::Index j=0;j<s.eigenvalues.size();++j){
     const double l=s.eigenvalues[j];
     modes<<j<<','<<l<<','<<std::sqrt(std::max(0.0,l))<<','
       <<(l<=s.zero_cutoff?"near_zero":"positive_candidate")<<','
       <<s.residuals[static_cast<std::size_t>(j)]<<'\n';
   }
   modes.close();
   const auto tm=rectangular_tm_reference(rect,4,4);
   std::ofstream bench(out/"tm_reference.csv");
   if(!bench)throw std::runtime_error("Cannot create TM reference CSV");
   bench<<std::setprecision(17)
       <<"m,n,lambda_analytic,kc_analytic,nearest_index,lambda_candidate,relative_error\n";
   std::vector<bool> used(static_cast<std::size_t>(s.eigenvalues.size()),false);
   for(const auto& t:tm){
     Eigen::Index best=-1;
     double distance=std::numeric_limits<double>::infinity();
     for(Eigen::Index j=0;j<s.eigenvalues.size();++j){
       if(used[static_cast<std::size_t>(j)]||
          s.eigenvalues[j]<=s.zero_cutoff)continue;
       const double d=std::abs(s.eigenvalues[j]-t.lambda);
       if(d<distance){distance=d;best=j;}
     }
     bench<<t.m<<','<<t.n<<','<<t.lambda<<','<<t.kc<<',';
     if(best>=0){
       used[static_cast<std::size_t>(best)]=true;
       bench<<best<<','<<s.eigenvalues[best]<<','<<distance/t.lambda;
     }else bench<<",,";
     bench<<'\n';
   }
   std::ofstream diag(out/"mode_diagnostics.csv");
   if(!diag)throw std::runtime_error("Cannot create diagnostic CSV");
   diag<<std::setprecision(17)
       <<"index,lambda,residual,divergence_ratio,normal_ratio,curl_norm_l2\n";
   int count=0;
   for(Eigen::Index j=0;j<s.eigenvalues.size()&&count<8;++j){
     const double l=s.eigenvalues[j];
     if(l<=s.zero_cutoff)continue;
     const auto d=modal_diagnostics(mls,quad,s.eigenvectors.col(j),l);
     diag<<j<<','<<l<<','
       <<s.residuals[static_cast<std::size_t>(j)]<<','
       <<d.divergence_ratio<<','<<d.normal_ratio<<','<<d.curl_l2<<'\n';
     if(count<4)write_field(out/("field_candidate_"+std::to_string(count+1)+".csv"),
                           mls,s.eigenvectors.col(j),rect,nx-1,ny-1);
     ++count;
   }
   diag.close();
   std::ofstream meta(out/"run_metadata.txt");
   meta<<std::setprecision(17)
       <<"model=cartesian_EFG_M3_experimental\n"
       <<"boundary=natural_curlH_zero\n"
       <<"warning=positive_candidates_are_not_validated_physical_modes\n"
       <<"nx="<<nx<<"\nny="<<ny<<"\ndegree="<<degree
       <<"\nsupport_factor="<<support
       <<"\nsupport_radius="<<radius
       <<"\nquadrature_order="<<quad.gauss_order
       <<"\nquadrature_points="<<a.quadrature_points
       <<"\nzero_cutoff="<<s.zero_cutoff
       <<"\nmass_condition="<<s.mass_condition
       <<"\nmax_mls_condition="<<a.max_moment_condition
       <<"\nnumerical_kernel="<<s.numerical_kernel_dimension<<'\n';
   SpectralExportOptions export_options;
   export_options.nx=nx;export_options.ny=ny;
   export_options.support_factor=support;
   export_options.origin="cartesian_spectrum_M3";
   export_spectrum(out,mls,quad,a,s,export_options);
   std::cout<<"M3 prototype: dofs="<<s.eigenvalues.size()
       <<" near_zero="<<s.numerical_kernel_dimension
       <<" mass_condition="<<s.mass_condition
       <<" output="<<out.string()<<'\n'
       <<"TM11 reference: "<<tm[0].lambda<<'\n'
       <<"CAUTION: physical modal convergence not yet established.\n";
 }catch(const std::exception& e){
   std::cerr<<"ERROR: "<<e.what()<<'\n';
   return 1;
 }
}
