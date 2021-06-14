#include<algorithm>
#include <math.h>

#include "tdrstyle.C"
#include "CMS_lumi.C"

class Likelihood
{
public:
    Float_t poi;
    Float_t deltaNLL;
    Likelihood(Float_t poi, Float_t deltaNLL): poi{poi}, deltaNLL{deltaNLL} {}

    bool operator<(Likelihood l) const
    {
        return this->poi < l.poi;
    }
};

double fa3_ggH_from_alpha_ggH(const double& alpha_ggH)
{

    //auto fa3 = pow(sin(alpha_ggH* M_PI / 180),2);
    auto fa3 = pow(sin(fabs(alpha_ggH)* M_PI / 180),2);

    if (alpha_ggH < 0.) fa3 *= -1;
    cout << alpha_ggH << "  " << fa3 << endl;
    return fa3;
}

double alpha_Hff_from_alpha_ggH(const double& alpha_ggH)
{
    //auto alpha = atan(tan(alpha_ggH* M_PI / 180)/2.38) * 180 / M_PI;
    auto alpha = atan(tan(alpha_ggH* M_PI / 180)/sqrt(2.38)) * 180 / M_PI;

    return alpha;
}

double alpha_ggH(const double& fa3_ggH)
{
    auto a1_ggH = sqrt(1.-fabs(fa3_ggH));
    auto a3_ggH = sqrt(fabs(fa3_ggH));

    auto alpha = atan(a3_ggH / a1_ggH) * 180 / M_PI;

    if (fa3_ggH < 0.) alpha *= -1;

    return alpha;
}

double alpha_Hff(const double& fa3_ggH)
{
    auto f = 1. / (1. + 2.38 * (1./fabs(fa3_ggH) -1.));

    auto alpha = asin(sqrt(f)) * 180 / M_PI;

    if (fa3_ggH < 0.) alpha *= -1;

    return alpha;
}

double reweight_fL1(const double& fL1)
{
    auto reweighted = 0.;
    auto sigma_r = 1.99;
    auto r = 1./(1.-2*0.23119); // r = 1/(cos(theta_w)^2 - sin(theta_w)^2) = 1/(1- 2*sin(theta_w)^2)
    if (fL1 != 0) reweighted = 1./(1.+(1./fabs(fL1)-1.)*(1.+pow(r,2)*sigma_r)/(1.+sigma_r))*fL1/fabs(fL1);

    return reweighted;
}

double reweight_fa2(const double& fa2)
{
    auto reweighted = 0.;
    auto r = 1.-0.23119; 
    auto sigma_r = 3.06;
    if (fa2 != 0) reweighted = 1./(1.+(1./fabs(fa2)-1.)*(1.+pow(r,2)*sigma_r)/(1.+sigma_r))*fa2/fabs(fa2);  

    return reweighted;
}

double reweight_fa3(const double& fa3)
{
    auto reweighted = 0.;
    auto r = 1.-0.23119;
    auto sigma_r = 3.15;
    if (fa3 != 0) reweighted = 1./(1.+(1./fabs(fa3)-1.)*(1.+pow(r,2)*sigma_r)/(1.+sigma_r))*fa3/fabs(fa3);

    return reweighted;
}

std::vector<Likelihood> read_scans(const bool& is_exp, const TString& coupling, const TString& channel, const TString& year)
{
    TString real_coupling = coupling;

    if (coupling == "alpha_ggH" || coupling == "alpha^Hff") real_coupling = "fa3_ggH";
    if (coupling == "alpha_ggH_ic" || coupling == "alpha^Hff_ic" || "fa3_ggH_ic") real_coupling = "alpha";

    if (coupling == "reweighted_fL1") real_coupling = "fL1";
    if (coupling == "reweighted_fa2") real_coupling = "fa2";
    if (coupling == "reweighted_fa3") real_coupling = "fa3";

    TString type = (is_exp) ? "exp" : "obs";

    TString input_file_location = "uscms_input/jan08/" + type + "/higgsCombine_scan_" + type + "_" + real_coupling + "_" + channel + "_" + year + ".MultiDimFit.mH125.root";
    if (coupling == "alpha_ggH_ic" || coupling == "alpha^Hff_ic" || "fa3_ggH_ic") {
      if(is_exp) input_file_location = "../output/paper_181120/cmb/125/higgsCombine.alpha.floatVBF.v2.exp.MultiDimFit.mH125.root";
      else input_file_location = "../output/paper_181120/cmb/125/higgsCombine.alpha.floatVBF.v2.obs.MultiDimFit.mH125.root";
   }

    auto input_file = std::make_unique<TFile>(input_file_location, "READ");
    std::unique_ptr<TTree> limit{static_cast<TTree*>(input_file->Get("limit"))};

    Float_t poi = 0.;
    Float_t deltaNLL = 0.;

    const auto poi_name = (real_coupling == "fa3_ggH" || real_coupling == "alpha") ? real_coupling : "CMS_zz4l_fai1";

    limit->SetBranchAddress(poi_name, &poi);
    limit->SetBranchAddress("deltaNLL", &deltaNLL);

    std::vector<Likelihood> scans;

    const auto& entries = limit->GetEntries();
    std::cout << "Read " << entries << " entries" << std::endl;

    auto count_best_fit = 0;

    for (auto entry = 0; entry < entries; ++entry)
    {
        limit->GetEntry(entry);

        if (deltaNLL == 0.)
        {
            if (count_best_fit == 0) count_best_fit++;
            else continue;
        }

        if (coupling == "alpha_ggH") poi = alpha_ggH(poi);
        else if (coupling == "alpha^Hff") poi = alpha_Hff(poi);
        else if (coupling == "reweighted_fL1") poi = reweight_fL1(poi);
        else if (coupling == "reweighted_fa2") poi = reweight_fa2(poi);
        else if (coupling == "reweighted_fa3") poi = reweight_fa3(poi);
        else if (coupling == "alpha^Hff_ic") poi = alpha_Hff_from_alpha_ggH(poi);
        else if (coupling == "fa3_ggH_ic") poi = fa3_ggH_from_alpha_ggH(poi);
        scans.push_back(Likelihood(poi, deltaNLL));
    }

    return scans;
}

TGraph* make_graph(const std::vector<Likelihood>& scans)
{
    const auto entries = scans.size();

    auto graph = new TGraph(0);
    auto poi = 0.;
    auto nll = 0.;

    for (auto entry = 0; entry < entries; ++entry)
    {
        poi = scans[entry].poi;
        nll = scans[entry].deltaNLL;

        graph->SetPoint(entry, poi, 2*nll);
    }

    return graph;
}

TString get_title(const TString& coupling)
{
    TString title = "";

    if (coupling == "fL1") title = "f_{#Lambda1}";
    else if (coupling == "fL1Zg") title = "f_{#Lambda1}^{Z#gamma}";
    else if (coupling == "fa2") title = "f_{a2}";
    else if (coupling == "fa3") title = "f_{a3}";
    else if (coupling == "alpha_ggH" || coupling == "alpha_ggH_ic") title = "#alpha_{ggH}";
    else if (coupling == "alpha^Hff" || coupling == "alpha^Hff_ic") title = "#alpha^{Hff}";
    else title = "f_{a3}^{ggH}";

    return title;
}

double get_min_poi(const TString& coupling)
{
    double min = 0.;

    if (coupling == "fa3_ggH" || coupling == "fa3_ggH_ic") min = -1.;
    else if(coupling == "alpha_ggH" || coupling == "alpha^Hff" || coupling == "alpha_ggH_ic" || coupling == "alpha^Hff_ic") min = -90.;
    else if(coupling == "reweighted_fL1") min = -0.04;
    else min = -0.1;

    return min; 
}

double get_max_poi(const TString& coupling)
{
    double max = 0.;

    if (coupling == "fa3_ggH" || coupling == "fa3_ggH_ic") max = 1.;
    else if(coupling == "alpha_ggH" || coupling == "alpha^Hff" || coupling == "alpha_ggH_ic" || coupling == "alpha^Hff_ic") max = 90.;
    else if(coupling == "reweighted_fL1") max = 0.04;
    else max = 0.1;

    return max; 
}

double get_max_2nll(const TString& coupling)
{
    double max = 0.;

    if (coupling == "fa3_ggH" || coupling == "alpha_ggH" || coupling == "alpha^Hff" || coupling == "alpha_ggH_ic" || coupling == "alpha^Hff_ic" || coupling == "fa3_ggH_ic") max = 10;
    else max = 45.;

    return max;
}

TString plot_name(const TString& coupling)
{
    TString plot_name = "";
    if (coupling == "fa3_ggH")
    {
        plot_name = "fa3ggH";
    }
    else if (coupling == "alpha_ggH")
    {
        plot_name = "alpha_ggH";
    }
    else if (coupling == "alpha^Hff")
    {
        plot_name = "alpha_Hff";
    }
    if (coupling == "fa3_ggH_ic")
    {   
        plot_name = "fa3_ggH_CB";
    }
    else if (coupling == "alpha_ggH_ic")
    {   
        plot_name = "alpha_ggH_CB";
    }
    else if (coupling == "alpha^Hff_ic")
    {   
        plot_name = "alpha_Hff_CB";
    }
    else
    {
        plot_name = coupling + "VBF";
    }
    return plot_name;
}

void plot_scan(const TString& coupling)
{
    auto CL68 = new TF1("CL68", "1", -1000000, 1000000);
    CL68->SetLineColor(1);
    CL68->SetLineStyle(9);

    auto CL95 = new TF1("CL95", "3.84", -1000000, 1000000);
    CL95->SetLineColor(1);
    CL95->SetLineStyle(9);

    auto scan_exp = read_scans(true, coupling, "emetmttt", "years");
    auto scan_obs = read_scans(false, coupling, "emetmttt", "years");

    sort(scan_exp.begin(), scan_exp.end());
    sort(scan_obs.begin(), scan_obs.end());

    std::cout << "Expected" << std::endl;
    for (auto scan : scan_exp)
    {
        std::cout << scan.poi << ": " << 2 * scan.deltaNLL << std::endl;
    }

    std::cout << "Observed" << std::endl;
    for (auto scan : scan_obs)
    {
        std::cout << scan.poi << ": " << 2 * scan.deltaNLL << std::endl;
    }

    auto curve_exp = make_graph(scan_exp);
    curve_exp->SetName("curve_exp");

    auto curve_obs = make_graph(scan_obs);
    curve_obs->SetName("curve_obs");

    curve_exp->SetLineColor(1);
    curve_obs->SetLineColor(1);

    curve_exp->SetLineWidth(3);
    curve_obs->SetLineWidth(3);

    curve_exp->SetLineStyle(2);

    auto coupling_title = get_title(coupling);

    auto min_poi = get_min_poi(coupling);
    auto max_poi = get_max_poi(coupling);
    auto min_2nll = -0.01;
    auto max_2nll = get_max_2nll(coupling);

    TH2D* frame = new TH2D("frame", "", 10, min_poi, max_poi, 20, min_2nll, max_2nll);
    frame->GetXaxis()->SetTitle(coupling_title);
    frame->GetXaxis()->CenterTitle();
    frame->GetXaxis()->SetTitleFont(43);
    frame->GetXaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetTitleSize(30);
    frame->GetXaxis()->SetLabelSize(0.03);

    frame->GetYaxis()->SetTitle("#minus 2 #Delta ln L");
    frame->GetYaxis()->CenterTitle();
    frame->GetYaxis()->SetTitleFont(43);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetYaxis()->SetTitleSize(30);
    frame->GetYaxis()->SetLabelSize(0.03);
    frame->SetStats(kFALSE);



    lumi_13TeV = "138 fb^{-1}";
    writeExtraText = (coupling == "alpha_ggH" || coupling == "alpha_ggH_ic") ? true : false;
    extraText = "       Supplementary";
    drawLogo = false;
    setTDRStyle();

    lumiTextSize = .4;
    cmsTextSize = .5;

    int W = 1200;
    int H = 1200;
    int H_ref = 1200;
    int W_ref = 1200;

    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;

    TCanvas* canvas = new TCanvas("canvas","",800,800);  
    canvas->SetFillColor(0);
    canvas->SetLeftMargin(L/W);
    canvas->SetRightMargin(R/W);
    canvas->SetTopMargin(T/H);
    canvas->SetBottomMargin(B/H);
    
    frame->Draw();


    if(coupling == "alpha^Hff_ic" || coupling == "alpha_ggH_ic" || coupling == "fa3_ggH_ic") {
      TSpline3 *spline_obs = new TSpline3("spline3_obs", curve_obs);
      TSpline3 *spline_exp = new TSpline3("spline3_exp", curve_exp);

      spline_exp->SetLineColor(1);
      spline_obs->SetLineColor(1);

      spline_exp->SetLineWidth(3);
      spline_obs->SetLineWidth(3);

      spline_exp->SetLineStyle(2);

      spline_obs->Draw("SAME");
      spline_exp->Draw("SAME");
    } else {
      curve_obs->Draw("CSAME");
      curve_exp->Draw("CSAME");
    }

    CL68->Draw("SAME");
    CL95->Draw("SAME");

    TLegend* legend = new TLegend(.15,.77,.44,.89);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(42);
    legend->AddEntry(curve_obs,"Observed","L");
    legend->AddEntry(curve_exp,"Expected","L");
    legend->SetBorderSize(0);
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.025);    

    latex.DrawLatex(0.23, 0.21,"68% CL");
    latex.DrawLatex(0.23, 0.44,"95% CL");

    CMS_lumi(canvas,4,0);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->GetFrame()->Draw();
    canvas->Print("plots/" + plot_name(coupling) + "_cmb.pdf");
}
