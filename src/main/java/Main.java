import com.vividsolutions.jts.geom.Geometry;

import edu.koritsas.slimemold.irrigation.IrrigationPhysarumPolycephalum;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumLagrarianCSPT;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumLangrarianCSP;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.*;
import org.geotools.graph.traverse.standard.DijkstraIterator;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) throws IOException {

       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
        //IrrigationNetwork network = new IrrigationNetwork("ParametrizedTests/HCSPT.shp", "ParametrizedTests/WCSPT.shp", "ParametrizedTests/PCSPT.shp");
        IrrigationNetwork network = new IrrigationNetwork("ParametrizedTests/Hbenchmark2.shp", "ParametrizedTests/Wbenchmark2.shp", "ParametrizedTests/Pbenchmark2.shp");
        //DirectedIrrigationNetwork network = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/HDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/WDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/PDMST.shp");
       Graph graph=null;
        try {
            graph =network.getBasicGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }


        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();




        Node source =sourceNodes.get(0);
        Node sink =sinkNodes.get(0);

      /*
         PhysarumPolycephalumLagrarianCSPT slime = new PhysarumPolycephalumLagrarianCSPT(graph,source,sinkNodes,1E-6,1E-6,200000,10000,1) {
             @Override
             public boolean pathViolatesConstraints(Graph graph) {

                 SimpleFeature sourcef = (SimpleFeature) source.getObject();
                 List<Edge> edgeList = new ArrayList<>(graph.getEdges());

                 double Ho= (double) sourcef.getAttribute("hdemand");

                 double dh = edgeList.stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
                     @Override
                     public double applyAsDouble(Edge edge) {


                         return getEdgeConstraintValue(edge);
                     }
                 }));



                 SimpleFeature sinkf = (SimpleFeature) sink.getObject();
                 double He= (double) sinkf.getAttribute("hdemand");

                 double DH= Ho-He;
                 System.out.println("Διαθέσιμο: "+DH+"Πραγματικό: "+dh+" Διαφορά: "+(DH-dh));
                            boolean violates=false;
                 if ((DH-dh)<0){
                     violates=true;

                 }



                return violates;
             }

             @Override
             public double getEdgeConstraintValue(Edge edge) {
                 SimpleFeature f = (SimpleFeature) edge.getObject();
                 //double dh = (double) f.getAttribute("Dh");

                 double D= (double) f.getAttribute("Diameter");
                 double Q= (double) f.getAttribute("Q");
                 Geometry g = (Geometry) f.getDefaultGeometry();
                 double L=g.getLength();

                 double dh= 0.00090940294*FastMath.pow(Q,1.78571428571)*L/FastMath.pow(D,4.78571428571);

                 return dh;
             }

             @Override
             public double getEdgeCost(Edge e) {
                 SimpleFeature f = (SimpleFeature) e.getObject();
                 Geometry g = (Geometry) f.getDefaultGeometry();
                 double L=g.getLength();
                 double cm = (double) f.getAttribute("Cost");
                 double cost=L*cm;
                 return cost;
             }
         };

           slime.execute();
        slime.showConductivityMap();
        slime.showFlowDiagram();

        GraphUtils.visualizeGraph(slime.getGraph());
        System.out.println(slime.pathViolatesConstraints(slime.getGraph()));
        System.out.println(slime.getSolutionCost());
     */


        IrrigationPhysarumPolycephalum slimeTree = new IrrigationPhysarumPolycephalum(graph,source,sinkNodes,1E-10,1E-10,300000,10000,1) {
            @Override
            public boolean pathViolatesConstraints(Graph graph) {

                DijkstraIterator.EdgeWeighter weigter = new DijkstraIterator.EdgeWeighter() {
                    @Override
                    public double getWeight(Edge e) {
                        return getEdgeConstraintValue(e);
                    }
                };


                DijkstraShortestPathFinder pf = new DijkstraShortestPathFinder(graph,sourceNode,weigter);
                pf.calculate();



                Node source=getSourceNode();


                SimpleFeature sourcef = (SimpleFeature) source.getObject();


                double Ho= (double) sourcef.getAttribute("hdemand");

                boolean violates =false;

                List<Edge> edgeList = new ArrayList<>(graph.getEdges());

               List<Node> hydr= network.getHydrants();
                for (Node n:hydr){
                    Path path = pf.getPath(n);

                    List<Edge> actual = new ArrayList<>();


                    for (Edge e:edgeList){
                        if (path.contains(e.getNodeA())&&path.contains(e.getNodeB())){
                            actual.add(e);
                        }
                    }


                    double dh = actual.stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
                        @Override
                        public double applyAsDouble(Edge edge) {


                            return getEdgeConstraintValue(edge);
                        }
                    }));



                    SimpleFeature sinkf = (SimpleFeature) n.getObject();
                    double He= (double) sinkf.getAttribute("hdemand");

                    double DH= Ho-He;
                    //System.out.println("Διαθέσιμο: "+DH+"Πραγματικό: "+dh+" Διαφορά: "+(DH-dh));

                    if (DH-dh<0){
                        violates=true;
                        break;
                    }
                }


                return violates;
            }

            @Override
            public double getEdgeConstraintValue(Edge edge) {
                SimpleFeature f = (SimpleFeature) edge.getObject();
                double dh = (double) f.getAttribute("Dh");

                double D= (double) f.getAttribute("Diameter");
                double Q= (double) f.getAttribute("Q");
                Geometry g = (Geometry) f.getDefaultGeometry();
                  double L=g.getLength();

               // double dh= 0.00090940294*FastMath.pow(Q,1.78571428571)*L/FastMath.pow(D,4.78571428571);

                return dh;
            }

            @Override
            public double getEdgeCost(Edge e) {
               SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry g = (Geometry) f.getDefaultGeometry();
                double L=g.getLength();
                double cm = (double) f.getAttribute("Cost");
                double cost=L*cm;
                return cost;

               // double Q= (double) f.getAttribute("Q");
                //double D = (double) f.getAttribute("Diameter");

                //double φ=FastMath.pow((FastMath.PI*D*D)/(4*9.77683),1/1.32559)*L*FastMath.pow(Q,0.4386);

               //double dh = (double) f.getAttribute("Dh");
                  //return 1/dh;
                //return φ;
            }
        };

        slimeTree.execute();
        slimeTree.showConductivityMap();
        slimeTree.showFlowDiagram();

      GraphUtils.visualizeGraph(slimeTree.getGraph());

        System.out.println(slimeTree.pathViolatesConstraints(slimeTree.getGraph()));

        System.out.println(slimeTree.getSolutionCost());

    }
}
