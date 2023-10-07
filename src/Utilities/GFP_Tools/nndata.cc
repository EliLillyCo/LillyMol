// Tools for dealing with NearNeighbours protos.

#include "nndata.h"

#include <stdio.h>

#include <iostream>

#include <google/protobuf/text_format.h>

#include "Foundational/iwstring/iwstring.h"

namespace gfp {

using std::cerr;

int
WriteNNData(const nnbr::NearNeighbours& proto, IWString_and_File_Descriptor& output) {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  if (!printer.PrintToString(proto, &buffer)) {
    cerr << "WriteNNData:cannot write '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  output << buffer;
  output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

}  // namespace gfp
