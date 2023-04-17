/* Fix to change Force computation in the Bohm GWP scheme to ensure energy conservation*/
#include "fix_Bohm_conserve.h"
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixConserve::FixConserve(LAMMPS *lmp, int narg, char **arg):
	Fix(lmp, narg, arg), id_temp(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal fix Conserve command");

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp";

}	

int FixConserve::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

void FixConserve::post_force()
{
	double **f = atom->f;
	double **beta = atom->beta_correction;

	int nlocal = atom->nlocal;
	if (igroup == atom->firstgroup) nlocal = atom->nfirst;

	for(int i = 0; i < nlocal; ++i){
		if (mask[i] & groupbit){
			

}