// SPDX-License-Identifier: MIT
// Al Maleh (2009), Chapter 6: fields built from gradients and z-cross-gradients.
// This driver emits all retained eigenpairs, including tiny/negative eigenvalues.
#include "al_maleh/rot_irrot.hpp"
#include <Eigen/Core>
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
struct Match { int nearest=-1,best=-1; double rho=0;};
void field_samples(const MLSEvaluator& mls,const RankReducedSpectrum& sp,
                   int selected,const std::filesystem::path& dest,
                   const std::string& component,double overlap) {
  if(selected<0)return;
  std::ofstream out(dest);
  if(!out)throw std::runtime_error("M5 cannot export H field");
  out<<std::setprecision(17)
     <<"sample_x,sample_y,hx,hy,selected_index,selected_lambda,shape_overlap\n";
  const Eigen::VectorXd v=sp.eigenvectors.col(selected);
  const auto n=static_cast<Eigen::Index>(mls.nodes().size());
  for(int j=0;j<=12;++j)for(int i=0;i<=24;++i) {
    const Point2 p{2.0*i/24.0,j/12.0};
    const auto f=evaluate_rot_irrot(mls.evaluate(p),v,n);
    double hx=f.hx,hy=f.hy;
    if(component=="gradient"){hx=f.irrot_x;hy=f.irrot_y;}
    else if(component=="rotational"){hx=f.rot_x;hy=f.rot_y;}
    if(!std::isfinite(hx)||!std::isfinite(hy))
      throw std::runtime_error("Nonfinite M5 shape field");
    out<<p.x<<','<<p.y<<','<<hx<<','<<hy<<','
       <<selected<<','<<sp.eigenvalues[selected]<<','<<overlap<<"\n";
  }
}
}
int main(int argc,char** argv) {
  try {
    if(argc<7||argc>10)
      throw std::invalid_argument(
        "Usage: rot_irrot_rect nx ny gauss_order support_factor subdivision output_dir "
        "[mass_rank_rel_cutoff=1e-10] [degree=1] [weight=wendland_c2|quartic_spline]");
    const int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]),gauss=std::stoi(argv[3]);
    const double factor=std::stod(argv[4]);
    const int sub=std::stoi(argv[5]);
    const std::filesystem::path outdir=argv[6];
    const double tolerance=argc>=8?std::stod(argv[7]):1e-10;
    const int degree=argc>=9?std::stoi(argv[8]):1;
    const MLSWeight weight=argc==10?parse_mls_weight(argv[9]):MLSWeight::WendlandC2;
    if(nx<3||ny<3||nx>13||ny>13||nx*ny>120||
       gauss<1||gauss>5||sub<1||sub>4||(degree!=1&&degree!=2)||
       !std::isfinite(factor)||factor<=0)
      throw std::invalid_argument("Invalid M5 experiment parameters");
    const Rectangle rect{0,2,0,1};
    const double radius=factor*std::max(2.0/(nx-1),1.0/(ny-1));
    const MLSEvaluator mls(make_regular_cloud(rect,nx,ny,radius),degree,1e12,weight);
    const BackgroundQuadrature quad{rect,(nx-1)*sub,(ny-1)*sub,gauss};
    const auto form=assemble_rot_irrot(mls,quad);
    const auto sol=solve_rot_irrot_rank_reduced(form,tolerance);
    std::filesystem::create_directories(outdir);
    const auto n=static_cast<Eigen::Index>(mls.nodes().size());
    {
      std::ofstream eig(outdir/"eigenvalues.csv");
      eig<<std::setprecision(17)
         <<"index,lambda,reduced_residual,full_residual,classification\n";
      for(Eigen::Index j=0;j<sol.eigenvalues.size();++j) {
        const double lambda=sol.eigenvalues[j];
        const char* category=lambda < -sol.eigenvalue_zero_cutoff ?
          "negative" : (std::abs(lambda)<=sol.eigenvalue_zero_cutoff?
            "near_zero":"positive");
        eig<<j<<','<<lambda<<','<<sol.reduced_residuals[static_cast<std::size_t>(j)]
           <<','<<sol.full_residuals[static_cast<std::size_t>(j)]
           <<','<<category<<'\n';
      }
    }
    {
      std::ofstream vec(outdir/"eigenvectors.csv");
      vec<<std::setprecision(17)<<"mode_index,dof,component,node,coefficient\n";
      for(Eigen::Index j=0;j<sol.eigenvectors.cols();++j)
        for(Eigen::Index i=0;i<sol.eigenvectors.rows();++i)
          vec<<j<<','<<i<<','<<(i<n?"gradient":"rotational")<<','
             <<i%n<<','<<sol.eigenvectors(i,j)<<'\n';
    }
    {
      std::ofstream matrices(outdir/"matrix_elements.csv");
      matrices<<std::setprecision(17)<<"row,col,K,M\n";
      for(Eigen::Index i=0;i<form.mass.rows();++i)
        for(Eigen::Index j=0;j<form.mass.cols();++j)
          matrices<<i<<','<<j<<','<<form.stiffness(i,j)<<','
                  <<form.mass(i,j)<<'\n';
    }
    // Original Chapter 6 Figs. 16–17: two vector shape functions from the
    // central scalar MLS basis function (5 rows × 9 columns of nodes).
    // Export our own field, never reuse the source author's image.
    if(nx==9 && ny==5) {
      const auto center=static_cast<std::size_t>((ny/2)*nx+nx/2);
      const char* names[]={"irrotational_shape_central.csv",
                           "rotational_shape_central.csv"};
      for(int kind=0;kind<2;++kind) {
        std::ofstream shapeout(outdir/names[kind]);
        if(!shapeout)throw std::runtime_error("Cannot export M5 central MLS shape");
        shapeout<<std::setprecision(17)
                <<"sample_x,sample_y,hx,hy,selected_index,selected_lambda,shape_overlap\n";
        for(int j=0;j<=12;++j)for(int i=0;i<=24;++i) {
          const Point2 p{2.0*i/24.0,j/12.0};
          const auto sh=mls.evaluate(p);
          double x=0,y=0;
          for(std::size_t k=0;k<sh.active_nodes.size();++k) {
            if(sh.active_nodes[k]!=center)continue;
            if(kind==0){x=sh.dx[k];y=sh.dy[k];}
            else{x=-sh.dy[k];y=sh.dx[k];}
          }
          shapeout<<p.x<<','<<p.y<<','<<x<<','<<y<<","
                  <<center<<",0,0\n";
        }
      }
    }
    const double pi=std::acos(-1.0);
    const double analytic[]={1.25*pi*pi,2*pi*pi};
    const int mm[]={1,2};
    std::vector<Point2> points;
    std::vector<double> weights;
    for_each_quadrature_point(quad,[&](Point2 p,double w){
      points.push_back(p);weights.push_back(w);
    });
    Eigen::MatrixXd curl=Eigen::MatrixXd::Zero(points.size(),sol.retained_rank);
    Eigen::VectorXd an1(points.size()),an2(points.size());
    for(std::size_t q=0;q<points.size();++q) {
      const auto sh=mls.evaluate(points[q]);
      for(std::size_t k=0;k<sh.active_nodes.size();++k) {
        const auto i=static_cast<Eigen::Index>(sh.active_nodes[k]);
        const double lap=sh.dxx[k]+sh.dyy[k];
        curl.row(q).noalias() += lap*sol.eigenvectors.row(n+i);
      }
      an1[static_cast<Eigen::Index>(q)]=
          std::sin(pi*points[q].x/2)*std::sin(pi*points[q].y);
      an2[static_cast<Eigen::Index>(q)]=
          std::sin(pi*points[q].x)*std::sin(pi*points[q].y);
    }
    std::vector<Match> matches(2);
    {
      std::ofstream out(outdir/"tm_matches.csv");
      out<<std::setprecision(17)
         <<"m,n,lambda_exact,nearest_index,nearest_lambda,nearest_signed_error_pct,"
         "best_curl_index,best_curl_lambda,best_curl_signed_error_pct,best_curl_rho,"
         "matched_by\n";
      for(int target=0;target<2;++target) {
        const Eigen::VectorXd& ref=target?an2:an1;
        double ref2=0;
        for(std::size_t k=0;k<weights.size();++k)
          ref2+=weights[k]*ref[static_cast<Eigen::Index>(k)]*
                            ref[static_cast<Eigen::Index>(k)];
        auto& m=matches[static_cast<std::size_t>(target)];
        double nearest=std::numeric_limits<double>::infinity();
        for(Eigen::Index j=0;j<sol.eigenvalues.size();++j) {
          const double lam=sol.eigenvalues[j];
          if(lam<=sol.eigenvalue_zero_cutoff)continue;
          const double gap=std::abs(lam-analytic[target]);
          if(gap<nearest){nearest=gap;m.nearest=static_cast<int>(j);}
          double curl2=0,mixed=0;
          for(std::size_t k=0;k<weights.size();++k) {
            const double c=curl(static_cast<Eigen::Index>(k),j);
            curl2+=weights[k]*c*c;
            mixed+=weights[k]*c*ref[static_cast<Eigen::Index>(k)];
          }
          const double overlap=(curl2>1e-15 && ref2>0) ?
            std::abs(mixed)/std::sqrt(curl2*ref2) : 0;
          if(overlap>m.rho){m.rho=overlap;m.best=static_cast<int>(j);}
        }
        if(m.nearest<0||m.best<0)
          throw std::runtime_error("M5 no positive TM candidate");
        auto error=[&](int idx){
          return 100.0*(sol.eigenvalues[idx]/analytic[target]-1);
        };
        out<<mm[target]<<",1,"<<analytic[target]<<','<<m.nearest<<','
           <<sol.eigenvalues[m.nearest]<<','<<error(m.nearest)<<','
           <<m.best<<','<<sol.eigenvalues[m.best]<<','
           <<error(m.best)<<','<<m.rho<<",spectral_and_curl_reported_separately\n";
        const auto label=target?"tm21":"tm11";
        field_samples(mls,sol,m.best,outdir/(std::string(label)+"_vector_field.csv"),
                      "total",m.rho);
        if(target==0) {
          field_samples(mls,sol,m.best,outdir/"tm11_irrot_field.csv",
                        "gradient",m.rho);
          field_samples(mls,sol,m.best,outdir/"tm11_rot_field.csv",
                        "rotational",m.rho);
        }
      }
    }
    {
      std::ofstream meta(outdir/"run_metadata.txt");
      meta<<std::setprecision(17)
          <<"method=Al_Maleh_2009_Ch6_Eqs6-1_to_6-6\n"
          <<"nx="<<nx<<"\nny="<<ny<<"\nnodes="<<n
          <<"\ndofs="<<2*n<<"\nquad_order="<<gauss
          <<"\nsubdivision="<<sub
          <<"\nmls_degree="<<degree
          <<"\nmls_weight="<<mls_weight_name(weight)
          <<"\nsupport_factor="<<factor
          <<"\nmass_rank_relative_cutoff="<<tolerance
          <<"\nmass_eigen_max_removed="<<sol.mass_eigen_max_removed
          <<"\nmass_kernel_stiffness_coupling="<<sol.mass_kernel_stiffness_coupling
          <<"\nmax_reduced_residual="<<sol.max_reduced_residual
          <<"\nmax_full_residual="<<sol.max_full_residual
          <<"\nintegration_points="<<form.quadrature_points
          <<"\nremoved_mass_nullity="<<sol.removed_mass_nullity
          <<"\nretained_rank="<<sol.retained_rank
          <<"\nmass_eigen_max="<<sol.mass_eigen_max
          <<"\nmass_min_retained="<<sol.mass_eigen_min_retained
          <<"\nmass_cutoff="<<sol.mass_rank_cutoff
          <<"\ncurl_zero_cutoff="<<sol.eigenvalue_zero_cutoff
          <<"\nnear_zero_eigenvalues="<<sol.numerical_curl_kernel
          <<"\nM_orthogonality_error="<<sol.M_orthogonality_error
          <<"\nmax_mls_condition="<<form.max_mls_condition
          <<"\nnegative_eigenvalues="
          <<(sol.eigenvalues.array()<-sol.eigenvalue_zero_cutoff).count()
          <<"\n";
    }
    std::cout<<std::setprecision(12)<<"M5 rot/irrot "<<ny<<"x"<<nx
             <<" gauss="<<gauss<<" rank="<<sol.retained_rank
             <<" mass_nullity="<<sol.removed_mass_nullity
             <<" zero="<<sol.numerical_curl_kernel
             <<" tm11_lambda="<<sol.eigenvalues[matches[0].best]
             <<" tm21_lambda="<<sol.eigenvalues[matches[1].best]
             <<" output="<<outdir.string()<<"\n";
  }catch(const std::exception& e){
    std::cerr<<"M5 ERROR: "<<e.what()<<"\n";return 1;
  }
}
