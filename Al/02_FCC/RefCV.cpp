/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/AtomValuePack.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>
#include <cfloat>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR REFCV
/*

*/
//+ENDPLUMEDOC


class RefCV : public MultiColvarBase {
private:
  double rcut2_, sigma_, sigmaSqr_;
  std::vector<std::vector<Vector>> Templates_;
  double lambda_;
public:
  static void registerKeywords( Keywords& keys );
  explicit RefCV(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Returns the number of coordinates of the field
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(RefCV,"REFCV")

void RefCV::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","SIGMA","0.1","Broadening parameter");
  keys.add("compulsory","CRYSTAL_STRUCTURE","FCC","Targeted crystal structure");
  keys.add("optional","LATTICE_CONSTANTS","Lattice constants");
  keys.add("compulsory","LAMBDA","100","Lambda parameter");
  keys.add("numbered","VECTOR","Relative distance of atoms from central atom.");
  keys.add("numbered","VECTOR1_","Relative distance of atoms from central atom of template 1.");
  keys.add("numbered","VECTOR2_","Relative distance of atoms from central atom of template 2.");
  keys.add("numbered","VECTOR3_","Relative distance of atoms from central atom of template 3.");
  keys.add("numbered","VECTOR4_","Relative distance of atoms from central atom of template 4.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
}

RefCV::RefCV(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  std::vector<double> lattice_constants;
  parseVector("LATTICE_CONSTANTS", lattice_constants);
  std::string crystal_structure;
  parse("CRYSTAL_STRUCTURE", crystal_structure);
  // find crystal structure
  double max_dist_ref_vector=0;
  if (crystal_structure == "FCC") {
    if (lattice_constants.size() != 1) error("Number of LATTICE_CONSTANTS arguments must be one for FCC");
    Templates_.resize(1);
    Templates_[0].resize(12);
    Templates_[0][0]  = Vector(+0.5,+0.5,+0.0)*lattice_constants[0];
    Templates_[0][1]  = Vector(-0.5,-0.5,+0.0)*lattice_constants[0];
    Templates_[0][2]  = Vector(+0.5,-0.5,+0.0)*lattice_constants[0];
    Templates_[0][3]  = Vector(-0.5,+0.5,+0.0)*lattice_constants[0];
    Templates_[0][4]  = Vector(+0.5,+0.0,+0.5)*lattice_constants[0];
    Templates_[0][5]  = Vector(-0.5,+0.0,-0.5)*lattice_constants[0];
    Templates_[0][6]  = Vector(-0.5,+0.0,+0.5)*lattice_constants[0];
    Templates_[0][7]  = Vector(+0.5,+0.0,-0.5)*lattice_constants[0];
    Templates_[0][8]  = Vector(+0.0,+0.5,+0.5)*lattice_constants[0];
    Templates_[0][9]  = Vector(+0.0,-0.5,-0.5)*lattice_constants[0];
    Templates_[0][10] = Vector(+0.0,-0.5,+0.5)*lattice_constants[0];
    Templates_[0][11] = Vector(+0.0,+0.5,-0.5)*lattice_constants[0];
    max_dist_ref_vector = std::sqrt(2)*lattice_constants[0]/2.;
  } else if (crystal_structure == "BCC") {
    if (lattice_constants.size() != 1) error("Number of LATTICE_CONSTANTS arguments must be one for BCC");
    Templates_.resize(1);
    Templates_[0].resize(14);
    Templates_[0][0]  = Vector(+0.5,+0.5,+0.5)*lattice_constants[0];
    Templates_[0][1]  = Vector(-0.5,-0.5,-0.5)*lattice_constants[0];
    Templates_[0][2]  = Vector(-0.5,+0.5,+0.5)*lattice_constants[0];
    Templates_[0][3]  = Vector(+0.5,-0.5,+0.5)*lattice_constants[0];
    Templates_[0][4]  = Vector(+0.5,+0.5,-0.5)*lattice_constants[0];
    Templates_[0][5]  = Vector(-0.5,-0.5,+0.5)*lattice_constants[0];
    Templates_[0][6]  = Vector(+0.5,-0.5,-0.5)*lattice_constants[0];
    Templates_[0][7]  = Vector(-0.5,+0.5,-0.5)*lattice_constants[0];
    Templates_[0][8]  = Vector(+1.0,+0.0,+0.0)*lattice_constants[0];
    Templates_[0][9]  = Vector(+0.0,+1.0,+0.0)*lattice_constants[0];
    Templates_[0][10] = Vector(+0.0,+0.0,+1.0)*lattice_constants[0];
    Templates_[0][11] = Vector(-1.0,+0.0,+0.0)*lattice_constants[0];
    Templates_[0][12] = Vector(+0.0,-1.0,+0.0)*lattice_constants[0];
    Templates_[0][13] = Vector(+0.0,+0.0,-1.0)*lattice_constants[0];
    max_dist_ref_vector = lattice_constants[0];
  } else if (crystal_structure == "CUSTOM") {
    max_dist_ref_vector = 0.0;
    // Case with one template 
    Templates_.resize(1);
    for(unsigned int i=1;; i++) {
      std::vector<double> tmp_vector(3);
      if(!parseNumberedVector("VECTOR",i,tmp_vector) ) {break;}
      Vector tmp_vector2; tmp_vector2[0]=tmp_vector[0]; tmp_vector2[1]=tmp_vector[1]; tmp_vector2[2]=tmp_vector[2];
      Templates_[0].push_back(tmp_vector2);
      double norm=tmp_vector2.modulo();
      if (norm>max_dist_ref_vector) max_dist_ref_vector=norm;
    }
    if (Templates_[0].size()==0) {
      // Case with multiple templates
      Templates_.clear();
      // Loop over templates
      for(unsigned int i=1;i<5; i++) {
        std::vector<Vector> tmp_template;
        std::string text = "VECTOR";
        text += std::to_string(i);
        text += "_";
        // Loop over template vectors
        for(unsigned int j=1;; j++) {
          std::vector<double> tmp_vector(3);
          if(!parseNumberedVector(text,j,tmp_vector) ) {break;}
          Vector tmp_vector2; tmp_vector2[0]=tmp_vector[0]; tmp_vector2[1]=tmp_vector[1]; tmp_vector2[2]=tmp_vector[2];
          tmp_template.push_back(tmp_vector2);
          double norm=tmp_vector2.modulo();
          if (norm>max_dist_ref_vector) max_dist_ref_vector=norm;
        }
        if (tmp_template.size()>0) Templates_.push_back(tmp_template);
      }
    }
    log.printf("  number of templates %d\n",Templates_.size() );
    log.printf("  number of vectors per template %d\n",Templates_[0].size() );
  } else {
    error("CRYSTAL_STRUCTURE=" + crystal_structure + " does not match any structures in the database");
  }
  log.printf("  targeting the %s crystal structure",crystal_structure.c_str());
  if (lattice_constants.size()>0) log.printf(" with lattice constants %f\n",lattice_constants[0]);
  else log.printf("\n");

  parse("SIGMA", sigma_);
  log.printf("  representing local density as a sum of Gaussians with standard deviation %f\n",sigma_);
  sigmaSqr_=sigma_*sigma_;
  log.printf("  maximum distance in the template is %f\n",max_dist_ref_vector);

  lambda_=100;
  parse("LAMBDA", lambda_);
  if (Templates_.size()>1) log.printf("  using a lambda value of %f\n",lambda_);

  // Set the link cell cutoff
  double rcut = max_dist_ref_vector + 3*sigma_;
  setLinkCellCutoff( rcut );
  rcut2_ = rcut * rcut;

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();
}

double RefCV::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  if (Templates_.size()==1) {
    // One template case
    for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
      Vector& distance=myatoms.getPosition(i);
      double d2;
      if ( (d2=distance[0]*distance[0])<rcut2_ &&
           (d2+=distance[1]*distance[1])<rcut2_ &&
           (d2+=distance[2]*distance[2])<rcut2_ &&
           d2>epsilon ) {
        // Iterate over atoms in the template
        for(unsigned k=0; k<Templates_[0].size(); ++k) {
          Vector distanceFromRef=distance-Templates_[0][k];
          double value = std::exp(-distanceFromRef.modulo2()/(4*sigmaSqr_) )/Templates_[0].size() ;
          // CAREFUL! Off-diagonal virial is incorrect. Do not perform NPT simulations with flexible box angles.
          accumulateSymmetryFunction( 1, i, value, (value/(2*sigmaSqr_))*(-distance+Templates_[0][k]) , (value/(2*sigmaSqr_))*Tensor(distance-Templates_[0][k],distance) , myatoms );
        }
      }
    }
    return myatoms.getValue(1);
  } else {
    // More than one template case
    std::vector<double> values(Templates_.size()); //value for each template
    // First time calculate sums
    for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
      Vector& distance=myatoms.getPosition(i);
      double d2;
      if ( (d2=distance[0]*distance[0])<rcut2_ &&
           (d2+=distance[1]*distance[1])<rcut2_ &&
           (d2+=distance[2]*distance[2])<rcut2_ &&
           d2>epsilon ) {
        // Iterate over templates
        for(unsigned j=0; j<Templates_.size(); ++j) {
          // Iterate over atoms in the template
          for(unsigned k=0; k<Templates_[j].size(); ++k) {
            Vector distanceFromRef=distance-Templates_[j][k];
            values[j] += std::exp(-distanceFromRef.modulo2()/(4*sigmaSqr_) )/Templates_[j].size() ;
          }
        }
      }
    }
    double sum=0;
    for(unsigned j=0; j<Templates_.size(); ++j) {
       values[j] = std::exp(lambda_*values[j]);
       sum += values[j];
    }
    // Second time find derivatives
    for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
      Vector& distance=myatoms.getPosition(i);
      double d2;
      if ( (d2=distance[0]*distance[0])<rcut2_ &&
           (d2+=distance[1]*distance[1])<rcut2_ &&
           (d2+=distance[2]*distance[2])<rcut2_ &&
           d2>epsilon ) {
        // Iterate over templates
        for(unsigned j=0; j<Templates_.size(); ++j) {
          // Iterate over atoms in the template
          for(unsigned k=0; k<Templates_[j].size(); ++k) {
            Vector distanceFromRef=distance-Templates_[j][k];
            double value = std::exp(-distanceFromRef.modulo2()/(4*sigmaSqr_) )/Templates_[j].size() ;
            accumulateSymmetryFunction( 1, i, value, -(values[j]/sum)*(value/(2*sigmaSqr_))*distanceFromRef  , (values[j]/sum)*(value/(2*sigmaSqr_))*Tensor(distanceFromRef,distance) , myatoms );
          }
        }
      }
    }
    return std::log(sum)/lambda_;
  }


}

}
}
