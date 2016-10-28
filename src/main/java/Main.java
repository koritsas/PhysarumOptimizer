import com.vividsolutions.jts.geom.Geometry;

import edu.koritsas.slimemold.mstree.PhysarumPolycephalumDirectedSPT;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumSPT;
import edu.koritsas.slimemold.shapefile.DirectedIrrigationNetwork;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.List;
import java.util.Random;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
        IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H10.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W10.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P10.shp");
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
       // Node sink =sinkNodes.get(random.nextInt(sinkNodes.size()));
/*
        Environment environment = new Environment(graph);
        FireAntColony<FireAnt> colony = new FireAntColony<FireAnt>(10) {
            @Override
            public FireAnt createFireAnt() {
                return new FireAnt(source) {
                    @Override
                    public boolean isEdgeValid(Edge e) {


                        return !(!getVisitedNodes().contains(e.getNodeA())&&!getVisitedNodes().contains(e.getNodeB()));
                    }

                    @Override
                    public List<Edge> getNeighbourhood(List<Edge> visitedEdges, List<Node> visitedNodes) {
                        List<Edge> neighbourhood = new ArrayList<>();
                        for (Node n:visitedNodes){
                            neighbourhood.addAll(n.getEdges());
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

                        return getSolution().getGraph().getNodes().containsAll(sinkNodes);
                    }

                    @Override
                    public double calculateSolutionCost(GraphBuilder solution) {
                       List<Edge> edges =new ArrayList<>(solution.getGraph().getEdges());
                        double cost =edges.stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
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
            }
        };

        AntSystemAlogrithm solver = new MaxMinAntSystemAlgorithm(environment,colony,300,8,0.5,0.5,0.5,0,15);
        solver.execute();

        GraphUtils.visualizeGraph(solver.getBestSolution().getGraph());
*/
        PhysarumPolycephalumSPT slimeMST = new PhysarumPolycephalumSPT(graph,source,sinkNodes,2,1.8,5000) {
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
