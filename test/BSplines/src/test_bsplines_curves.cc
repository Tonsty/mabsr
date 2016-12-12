#include <MABSR/BSplines.h>

#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkAxis.h>

extern int verbose = 0;

int main(int, char *[])
{
	// Create a table with some points in it
	vtkSmartPointer<vtkTable> table = 
		vtkSmartPointer<vtkTable>::New();

	vtkSmartPointer<vtkFloatArray> arrX1 = 
		vtkSmartPointer<vtkFloatArray>::New();
	arrX1->SetName("X Axis 1");
	table->AddColumn(arrX1);

	vtkSmartPointer<vtkFloatArray> arrB1 = 
		vtkSmartPointer<vtkFloatArray>::New();
	arrB1->SetName("B-Spline-order-1");
	table->AddColumn(arrB1);

	vtkSmartPointer<vtkFloatArray> arrX2 = 
		vtkSmartPointer<vtkFloatArray>::New();
	arrX2->SetName("X Axis 2");
	table->AddColumn(arrX2);

	vtkSmartPointer<vtkFloatArray> arrB2 = 
		vtkSmartPointer<vtkFloatArray>::New();
	arrB2->SetName("B-Spline-order-2");
	table->AddColumn(arrB2);

	vtkSmartPointer<vtkFloatArray> arrX3 = 
		vtkSmartPointer<vtkFloatArray>::New();
	arrX3->SetName("X Axis 3");
	table->AddColumn(arrX3);

	vtkSmartPointer<vtkFloatArray> arrB3 = 
		vtkSmartPointer<vtkFloatArray>::New();
	arrB3->SetName("B-Spline-order-3");
	table->AddColumn(arrB3);

	// Fill in the table with some example values
	
	float X_length = 2.0f;

	//int numKnots = 8; 
	//mabsr::Vector knots(numKnots);
	//float inc_knots = 1.0 / (numKnots - 1);
	//for (int i = 0; i < numKnots; i++) {
	//	knots(i) = i * inc_knots;
	//}
	//mabsr::BSplines bsplines(knots);

	//int numPoints = 400;
	//float inc = 1.0 / (numPoints-1);
	//table->SetNumberOfRows(numPoints);
	//for (int i = 0; i < numPoints; ++i)
	//{
	//	table->SetValue(i, 0, i * inc - 0.5);
	//	table->SetValue(i, 1, bsplines.B(0, 5, i * inc));
	//	table->SetValue(i, 2, i * inc - 0.5);
	//	table->SetValue(i, 3, bsplines.B(1, 5, i * inc));
	//	table->SetValue(i, 4, i * inc - 0.5);
	//	table->SetValue(i, 5, bsplines.B(2, 5, i * inc));
	//}

	int numq = 32;
	mabsr::Matrix q(numq, 2);
	float inc_controls = X_length / (numq - 1);
	mabsr::Vector random_Y = mabsr::Vector::Random(numq);
	for(int i = 0; i < numq; i++) {
		q(i, 0) = i * inc_controls - 1.0f;
		q(i, 1) = random_Y(i) * 0.2 + q(i,0) * ((i < numq/2) ? -1 : 1);
	}

	if(verbose>=1) std::cerr << "q = \n" << q << std::endl;

	mabsr::InterpolationBsplines inter_bsplines_1(q, 1);
	mabsr::InterpolationBsplines inter_bsplines_2(q, 4, true);
	//mabsr::InterpolationBsplines inter_bsplines_3(q, 3, true);
	mabsr::ApproximationBsplines approxi_bsplines(q, 4, 30, 0.0000001);

	int numPoints = 800;
	float inc = 0.999999f   / (numPoints-1);
	table->SetNumberOfRows(numPoints);
	for (int i = 0; i < numPoints; ++i)
	{
		mabsr::Point point1 = inter_bsplines_1.eval(i * inc + 0.0000005f );
		table->SetValue(i, 0, point1(0));
		table->SetValue(i, 1, point1(1));
		mabsr::Point point2 = inter_bsplines_2.eval(i * inc + 0.0000005f);
		table->SetValue(i, 2, point2(0));
		table->SetValue(i, 3, point2(1));
		//mabsr::Point point3 = inter_bsplines_3.eval(i * inc + 0.0005f);etValue(i, 4, point3(0));
		//table->SetValue(i, 5, point3(1));
	}
	float inc2 = 1.0f / (numPoints-1);
	for (int i = 0; i < numPoints; ++i)
	{
		mabsr::Point point3 = approxi_bsplines.eval(i * inc2 - (i == numPoints-1) * 0.00000005);
		table->SetValue(i, 4, point3(0));
		table->SetValue(i, 5, point3(1));
	}

	// Set up the view
	vtkSmartPointer<vtkContextView> view = 
		vtkSmartPointer<vtkContextView>::New();
	view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

	// Add multiple line plots, setting the colors etc
	vtkSmartPointer<vtkChartXY> chart = 
		vtkSmartPointer<vtkChartXY>::New();
	chart->GetAxis(vtkAxis::BOTTOM)->SetRange(-1.2, 1.2);
	chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::FIXED);
	chart->GetAxis(vtkAxis::LEFT)->SetRange(-1.2, 1.2);
	chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
	view->GetScene()->AddItem(chart);                                                                                                 
	vtkPlot *line = chart->AddPlot(vtkChart::LINE);
#if VTK_MAJOR_VERSION <= 5
	line->SetInput(table, 4, 5);
#else
	line->SetInputData(table, 0, 1);
#endif
	line->SetColor(0, 255, 0, 255);
	line->SetWidth(3.0);
	line = chart->AddPlot(vtkChart::LINE);
#if VTK_MAJOR_VERSION <= 5
	line->SetInput(table, 2, 3);
#else
	line->SetInputData(table, 2, 3);
#endif
	line->SetColor(255, 0, 0, 255);
	line->SetWidth(2.0);
	line = chart->AddPlot(vtkChart::POINTS);
#if VTK_MAJOR_VERSION <= 5
	line->SetInput(table, 0, 1);
#else
	line->SetInputData(table, 4, 5);
#endif
	line->SetColor(0, 0, 255, 255);
	line->SetWidth(1.0);

	// For dotted line, the line type can be from 2 to 5 for different dash/dot
	// patterns (see enum in vtkPen containing DASH_LINE, value 2):
#ifndef WIN32
	line->GetPen()->SetLineType(vtkPen::DASH_LINE);
#endif
	// (ifdef-ed out on Windows because DASH_LINE does not work on Windows
	//  machines with built-in Intel HD graphics card...)

	//view->GetRenderWindow()->SetMultiSamples(0);
	view->GetRenderWindow()->SetSize(600, 600);

	// Start interactor
	view->GetInteractor()->Initialize();
	view->GetInteractor()->Start();

	return EXIT_SUCCESS;
}