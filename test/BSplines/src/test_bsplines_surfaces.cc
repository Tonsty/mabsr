#include <MABSR/BSplines.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>

extern int verbose = 0;

int main(int, char *[])
{
	float X_length = 2.0f, Y_length = 2.0f;

	int numq_x = 32, numq_y = 32;
	mabsr::PointSet q(numq_x * numq_y, 3);
	float inc_controls_x = X_length / (numq_x - 1), inc_controls_y = Y_length / (numq_y - 1);
	mabsr::Vector random_Z = mabsr::Vector::Random(numq_x * numq_y);
	for(int i = 0; i < numq_y; i++) {
		for (int j = 0; j < numq_x; j++) {
			q(i * numq_y + j, 0) = i * inc_controls_y - 1.0f;
			q(i * numq_y + j, 1) = j * inc_controls_x - 1.0f;
			q(i * numq_y + j, 2) = random_Z(i * numq_y + j) * 0.05+ (q(i * numq_y + j, 0) * ((i < numq_y/2) ? -1 : 1) + q(i * numq_y + j, 1) * ((j < numq_x/2) ? -1 : 1)) * 0.5;
		}
	}

	mabsr::Vector bar_u(numq_x), bar_v(numq_y);
	for (int j = 0; j < numq_x; j++){
		bar_u(j) = j * 1.0/(numq_x-1);
	}
	for (int i = 0; i < numq_y; i++) {
		bar_v(i) = i * 1.0/(numq_y-1);
	}

	//mabsr::TensorBSplines tensorbsplines(numq_x, numq_y, q, bar_u, bar_v, 3, 3, 16, 16);

	mabsr::Vector bar_u_new(numq_x * numq_y), bar_v_new(numq_x * numq_y);
	for(int j = 0; j < numq_y; j++) {
		for(int i = 0; i < numq_x; i++) {
			bar_u_new(j * numq_x + i) = bar_u(i);
			bar_v_new(j * numq_x + i) = bar_v(j);
		}
	}

	mabsr::adapt_bar_t(bar_u_new);
	mabsr::adapt_bar_t(bar_v_new);

	if(verbose>=1) std::cout << "bar_u_new = \n" << bar_u_new << std::endl;
	if(verbose>=1) std::cout << "bar_v_new = \n" << bar_v_new << std::endl;

	mabsr::TensorBSplines tensorbsplines(numq_x * numq_y, q, bar_u_new, bar_v_new, 3, 3, 16, 16, 0.001);	

	// Setup points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkPoints> points1 =
		vtkSmartPointer<vtkPoints>::New();

	for (int i = 0; i < numq_x * numq_y; i++) {
		points->InsertNextPoint(q(i, 0), q(i, 1), q(i, 2));
	}

	int res_x = 100, res_y = 100;

	for (int i = 0; i < res_y; i++) {
		for(int j = 0; j < res_x; j++) {
			mabsr::Point p = tensorbsplines.eval( j * 1.0 / (res_x-0.9999), i * 1.0 /(res_y-0.9999));
			points1->InsertNextPoint(p(0), p(1), p(2));
		}
	}

	vtkSmartPointer<vtkCellArray> polygons =
		vtkSmartPointer<vtkCellArray>::New();

	vtkSmartPointer<vtkCellArray> polygons1 =
		vtkSmartPointer<vtkCellArray>::New();

	for (int i = 0; i < numq_y - 1; i++) {
		for (int j = 0; j < numq_x - 1; j++) {
			// Create the polygon
			vtkSmartPointer<vtkPolygon> polygon =
				vtkSmartPointer<vtkPolygon>::New();
			polygon->GetPointIds()->SetNumberOfIds(4); //make a quad
			polygon->GetPointIds()->SetId(0, i * numq_x + j);
			polygon->GetPointIds()->SetId(1, i * numq_x + j + 1);
			polygon->GetPointIds()->SetId(2, (i + 1) * numq_x + j + 1);
			polygon->GetPointIds()->SetId(3, (i + 1) * numq_x + j);

			// Add the polygon to a list of polygons
			polygons->InsertNextCell(polygon);
		}
	}

	for (int i = 0; i < res_y-1; i++) {
		for (int j = 0; j < res_x-1; j++) {
			// Create the polygon
			vtkSmartPointer<vtkPolygon> polygon =
				vtkSmartPointer<vtkPolygon>::New();
			polygon->GetPointIds()->SetNumberOfIds(4); //make a quad
			polygon->GetPointIds()->SetId(0, i * res_x + j);
			polygon->GetPointIds()->SetId(1, i * res_x + j + 1);
			polygon->GetPointIds()->SetId(2, (i + 1) * res_x + j + 1);
			polygon->GetPointIds()->SetId(3, (i + 1) * res_x + j);

			// Add the polygon to a list of polygons
			polygons1->InsertNextCell(polygon);
		}
	}

	// Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(points);
	polygonPolyData->SetPolys(polygons);

	vtkSmartPointer<vtkPolyData> polygonPolyData1 =
		vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData1->SetPoints(points1);
	polygonPolyData1->SetPolys(polygons1);

	// Create a mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(polygonPolyData);
#else
	mapper->SetInputData(polygonPolyData);
#endif

	vtkSmartPointer<vtkPolyDataMapper> mapper1 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper1->SetInput(polygonPolyData1);
#else
	mapper1->SetInputData(polygonPolyData1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	actor->GetProperty()->SetColor(0.0, 0.0, 1.0);
	actor->GetProperty()->SetRepresentationToWireframe();

	vtkSmartPointer<vtkActor> actor1 =
		vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper1);

	//actor1->GetProperty()->SetRepresentationToWireframe();

	// Visualize
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

	renderWindowInteractor->SetInteractorStyle( style );

	renderer->AddActor(actor);
	renderer->AddActor(actor1);
	renderer->SetBackground(.5,.3,.31); // Background color salmon
	renderWindow->SetSize(600, 600);

	vtkSmartPointer<vtkAxesActor> axes = 
		vtkSmartPointer<vtkAxesActor>::New();

	vtkSmartPointer<vtkOrientationMarkerWidget> widget = 
		vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
	widget->SetOrientationMarker( axes );
	widget->SetInteractor( renderWindowInteractor );
	widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
	widget->SetEnabled( 1 );
	widget->InteractiveOn();

	renderer->ResetCamera();
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}