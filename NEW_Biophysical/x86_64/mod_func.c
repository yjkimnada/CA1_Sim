#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _CaT_reg(void);
extern void _ca_reg(void);
extern void _cad_reg(void);
extern void _cagk_reg(void);
extern void _carF_reg(void);
extern void _exp2synNMDA_reg(void);
extern void _id_reg(void);
extern void _kad_reg(void);
extern void _kap_reg(void);
extern void _kca_reg(void);
extern void _kdr_reg(void);
extern void _km_reg(void);
extern void _kv_reg(void);
extern void _na_reg(void);
extern void _nadend_reg(void);
extern void _nax_reg(void);
extern void _vecevent_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," .//CaT.mod");
    fprintf(stderr," .//ca.mod");
    fprintf(stderr," .//cad.mod");
    fprintf(stderr," .//cagk.mod");
    fprintf(stderr," .//carF.mod");
    fprintf(stderr," .//exp2synNMDA.mod");
    fprintf(stderr," .//id.mod");
    fprintf(stderr," .//kad.mod");
    fprintf(stderr," .//kap.mod");
    fprintf(stderr," .//kca.mod");
    fprintf(stderr," .//kdr.mod");
    fprintf(stderr," .//km.mod");
    fprintf(stderr," .//kv.mod");
    fprintf(stderr," .//na.mod");
    fprintf(stderr," .//nadend.mod");
    fprintf(stderr," .//nax.mod");
    fprintf(stderr," .//vecevent.mod");
    fprintf(stderr, "\n");
  }
  _CaT_reg();
  _ca_reg();
  _cad_reg();
  _cagk_reg();
  _carF_reg();
  _exp2synNMDA_reg();
  _id_reg();
  _kad_reg();
  _kap_reg();
  _kca_reg();
  _kdr_reg();
  _km_reg();
  _kv_reg();
  _na_reg();
  _nadend_reg();
  _nax_reg();
  _vecevent_reg();
}
