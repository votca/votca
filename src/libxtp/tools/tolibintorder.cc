#include "tolibintorder.h"


namespace votca {
namespace xtp {

void ToLibintOrder::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
}

bool ToLibintOrder::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  
  XTP_LOG(Log::error, _log)
      << "\n Some MESSAGE \n";
  
  return true;
}

}  // namespace xtp
}  // namespace votca