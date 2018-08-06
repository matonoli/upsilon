//
// Pibero Djawotho <pibero@indiana.edu>
// Indiana University
// 3 June 2008
//

#include "StUpsilonCandidate.h"

ClassImp(StUpsilonCandidate);

StUpsilonCandidate::StUpsilonCandidate(const StUpsilonCandidate& ups)
{
  if (this == &ups) return;
  *this = ups;
}
