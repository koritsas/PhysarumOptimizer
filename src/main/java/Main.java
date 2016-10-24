import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.aco.components.ACOProblemSolver;
import edu.koritsas.aco.components.Environment;
import edu.koritsas.aco.components.FireAnt;
import edu.koritsas.aco.components.FireAntColony;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumDirectedMST;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumMST;
import edu.koritsas.slimemold.shapefile.DirectedIrrigationNetwork;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.graph.build.GraphBuilder;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
        IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P4.shp");
        //DirectedIrrigationNetwork network = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/HDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/WDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/PDMST.shp");
       Graph graph=null;
        try {
            graph =network.getBasicGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }


        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();
        //sinkNodes=sinkNodes.subList(7,8);

        Node source =sourceNodes.get(0);

        Random random = new Random();
        Node sink =sinkNodes.get(random.nextInt(sinkNodes.size()));


        FireAnt fireAnt =new FireAnt(source) {
            @Override
            public List<Edge> getNeighbourhood(List<Edge> visitedEdges, List<Node> visitedNodes) {
                List<Edge> neighbourhood = new ArrayList<>();
                synchronized (neighbourhood) {
                    for (Node n : visitedNodes) {
                        neighbourhood.addAll(n.getEdges());
                    }
                }
                return neighbourhood;
            }

            @Override
            public boolean violatesConstraints(Edge edge) {
                return false;
            }

            @Override
            public double getEdgeHeuristicValue(Edge edge) {
                SimpleFeature f = (SimpleFeature) edge.getObject();
                Geometry g = (Geometry) f.getDefaultGeometry();
                return 1/g.getLength();
            }

            @Override
            public boolean isSolutionCompleted() {

                return (getSolution().getGraph().getNodes().size()-1)==sinkNodes.size();
            }

            @Override
            public double calculateSolutionCost(GraphBuilder solution) {
                double cost= (double) solution.getGraph().getEdges().stream().filter(o -> o==null).collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
                    @Override
                    public double applyAsDouble(Edge edge) {
                        SimpleFeature f = (SimpleFeature) edge.getObject();
                        Geometry g = (Geometry) f.getDefaultGeometry();
                        return g.getLength();
                    }
                }));


                return cost;
            }
        };
        Environment environment = new Environment(graph);
        FireAntColony fireAntColony = new FireAntColony(10) {
            @Override
            public FireAnt createFireAnt() {
                return fireAnt;
            }
        };

      ACOProblemSolver  solver = new ACOProblemSolver(environment,fireAntColony,100,1,0.5,0.5,0.5);
         solver.execute();
        GraphUtils.visualizeGraph(solver.getBestSolution().getGraph());


        PhysarumPolycephalumMST slimeMST = new PhysarumPolycephalumMST(graph,source,sinkNodes,1,1.8,50000) {
            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry g = (Geometry) f.getDefaultGeometry();
                return g.getLength();
            }
        };

        slimeMST.execute();
        slimeMST.showFlowDiagram();
        slimeMST.showConductivityMap();
        slimeMST.showPressureMap();


        GraphUtils.visualizeGraph(graph);

        System.out.println(slimeMST.getSolutionCost());

    }
}
