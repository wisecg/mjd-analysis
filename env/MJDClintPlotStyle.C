// Clint's settings.  Adapted from the MJDTalkPlotStyle.C
//
{
  //WARNING: If histograms not made with this style
  //You will need gROOT->ForceStyle() in your macro.
  printf("Loading the Clint Plot Style\n");

  // Title, Stats, Date off by default
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);

  // Other text (e.g. legends, stats, histo title)
  kTextFont = 49; // 4 (Helvetica) + 9 (size in pixels)
  kTextSize = 30;
  gStyle->SetTextSize(kTextSize);
  gStyle->SetTextFont(kTextFont);
  gStyle->SetTextColor(1);

  // Axis titles
  kAxisTitleFont = 43; // 4 (Helvetica) + 3 (size in pixels)
  kAxisTitleSize = 30;
  kAxisTitleOffset = 1;
  gStyle->SetTitleSize(kAxisTitleSize, "XYZ");
  gStyle->SetTitleFont(kAxisTitleFont, "XYZ");
  gStyle->SetTitleXOffset(kAxisTitleOffset);
  gStyle->SetTitleYOffset(kAxisTitleOffset);

  // Axis labels
  kLabelFont = 43;
  kLabelSize = 20;	// 35 for blind people
  kLabelOffset = 0.006;
  kLabelSizeZ = 20;
  gStyle->SetLabelFont(kLabelFont, "XYZ");
  gStyle->SetLabelSize(kLabelSize, "XYZ");
  gStyle->SetLabelSize(kLabelSizeZ, "Z");
  gStyle->SetLabelOffset(kLabelOffset, "XY");
  gStyle->SetLabelOffset(kLabelOffset*0.5, "Z");


  // Fill solid by default
  kFillStyle=1001;

  // No little lines at the ends of error bars.
  gStyle->SetEndErrorSize(0);

  // Canvas width and height: 600x800
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(800);
  gStyle->SetCanvasBorderMode(0);

  // Pads and margins:
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);

  //Frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineWidth(1);

  // // Legend: borderless with white background.
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

  // // Colors: use rainbow scale
  gStyle->SetPalette(1);
}

