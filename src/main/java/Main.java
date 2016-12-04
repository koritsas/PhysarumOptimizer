import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumLagrarianCSPT;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.graph.traverse.standard.DijkstraIterator;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
        IrrigationNetwork network = new IrrigationNetwork("ParametrizedTests/Hbenchmark2.shp", "ParametrizedTests/Wbenchmark2.shp", "ParametrizedTests/Pbenchmark2.shp");
        //DirectedIrrigationNetwork network = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/HDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/WDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/PDMST.shp");
       Graph graph=null;
        try {
            graph =network.getBasicGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }

   /*     graph.getEdges().removeIf(new Predicate() {
            @Override
            public boolean test(Object o) {
                Edge e= (Edge) o;
                SimpleFeature f = (SimpleFeature) e.getObject();
                double v = (double) f.getAttribute("V");
                if (v>2){
                    return true;
                }else {
                    return false;
                }
            }
        });*/

        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();
        //sinkNodes=sinkNodes.subList(7,8);

        Node source =sourceNodes.get(0);
     /*
        PhysarumPolycephalumLangrarianCSP slime = new PhysarumPolycephalumLangrarianCSP(graph,source,sinkNodes.get(0),1E-12,1E-12,50000,100,0.1) {
            @Override
            public boolean pathViolatesConstraints(Graph graph) {
                Node source=getSourceNode();
                Node sink = getSinkNode();

                SimpleFeature sourcef = (SimpleFeature) source.getObject();
                SimpleFeature sinkf = (SimpleFeature) sink.getObject();

                double Ho= (double) sourcef.getAttribute("hdemand");
                double He= (double) sinkf.getAttribute("hdemand");

                double dh = (double) graph.getEdges().stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
                    @Override
                    public double applyAsDouble(Edge edge) {


                        return getEdgeConstraintValue(edge);
                    }
                }));
                boolean violates=false;

                if ((Ho-He)<dh){
                    violates=true;
                }

                return violates;
            }

            @Override
            public double getEdgeConstraintValue(Edge edge) {
               SimpleFeature f = (SimpleFeature) edge.getObject();
                double dh = (double) f.getAttribute("Dh");
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
        PhysarumPolycephalumLagrarianCSPT slimeTree = new PhysarumPolycephalumLagrarianCSPT(graph,source,sinkNodes,1E-12,1E-12,25000,100,0.1) {
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

                for (Node n:sinkNodes){
                    Path path = pf.getPath(n);
                    List<Edge> edges = path.getEdges();

                    double dh = (double) graph.getEdges().stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
                        @Override
                        public double applyAsDouble(Edge edge) {


                            return getEdgeConstraintValue(edge);
                        }
                    }));



                    SimpleFeature sinkf = (SimpleFeature) n.getObject();
                    double He= (double) sinkf.getAttribute("hdemand");

                    if ((Ho-He)<dh){
                        violates=true;
                    }
                }


                return violates;
            }

            @Override
            public double getEdgeConstraintValue(Edge edge) {
                SimpleFeature f = (SimpleFeature) edge.getObject();
                double dh = (double) f.getAttribute("Dh");
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

        slimeTree.execute();
        slimeTree.showConductivityMap();
        slimeTree.showFlowDiagram();

        GraphUtils.visualizeGraph(slimeTree.getGraph());

        System.out.println(slimeTree.pathViolatesConstraints(slimeTree.getGraph()));

        System.out.println(slimeTree.getSolutionCost());
    }
}
