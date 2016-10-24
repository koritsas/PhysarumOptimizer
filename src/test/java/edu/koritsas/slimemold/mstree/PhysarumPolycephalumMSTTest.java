package edu.koritsas.slimemold.mstree;

import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.aco.components.ACOProblemSolver;
import edu.koritsas.aco.components.Environment;
import edu.koritsas.aco.components.FireAnt;
import edu.koritsas.aco.components.FireAntColony;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.graph.build.GraphBuilder;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 24/10/2016.
 */
@RunWith(Parameterized.class)
public class PhysarumPolycephalumMSTTest {

    public PhysarumPolycephalumMSTTest(IrrigationNetwork network){
        this.network=network;
    }
    private final Random random=new Random();


    private IrrigationNetwork network;

    private final  double Io=2;
    private final double γ=1.8;
    private final int numberOfIterations=300;
    private PhysarumPolycephalumMST slime;
    private Graph graph;
    private DijkstraShortestPathFinder pf;
    private Node source;
    private Node sink;
    private List<Node> sinkNodes;


    private FireAnt fireAnt;
    private Environment environment;
    private FireAntColony fireAntColony;
    private ACOProblemSolver solver;

    @Before
    public void setUp() throws Exception {

        Graph graph=network.getBasicGraph();
        List<Node> v=network.getWaterSource();
        source=network.getWaterSource().get(0);;
        //sink=network.getHydrants().get(random.nextInt(network.getHydrants().size()));
        sinkNodes =network.getHydrants();
        slime = new PhysarumPolycephalumMST(graph,source,sinkNodes,Io,γ,numberOfIterations) {
            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry geometry = (Geometry) f.getDefaultGeometry();
                return geometry.getLength();
            }
        };

        fireAnt = new FireAnt(source) {
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

                return getVisitedNodes().containsAll(sinkNodes);
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
         environment = new Environment(graph);
        fireAntColony = new FireAntColony(10) {
            @Override
            public FireAnt createFireAnt() {
                return fireAnt;
            }
        };

       solver = new ACOProblemSolver(environment,fireAntColony,100,1,0.5,0.5,0.5);



    }

    @After
    public void tearDown() throws Exception {
        network=null;

        slime=null;
        source=null;
        sink=null;
        sinkNodes=null;
        fireAnt=null;
        fireAntColony=null;
        solver=null;
        environment=null;
    }
    @Parameterized.Parameters
    public static Collection graphs() throws IOException {
        IrrigationNetwork network1 = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P1.shp");

        IrrigationNetwork network2 = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H2.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W2.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P2.shp");

        IrrigationNetwork network3 = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H3.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W3.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P3.shp");

        IrrigationNetwork network4 = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P4.shp");
        ////Graph g1 = network1.getGraph();
        // Graph g2 = network2.getGraph();
        // Graph g3 = network3.getGraph();


        return Arrays.asList(new IrrigationNetwork[]{network1, network2, network3,network4});
    };

    @Test
    public void execute() throws Exception {
        solver.execute();
        GraphBuilder antGraphBuilder=solver.getBestSolution();

        slime.execute();
        slime.getGraph();

        Assert.assertTrue(solver.getBestSolutionCost()==slime.getSolutionCost());
        Assert.assertTrue(solver.getBestSolution().getGraph().getEdges().containsAll(slime.getGraph().getEdges()));
        Assert.assertTrue(slime.getGraph().getEdges().containsAll(solver.getBestSolution().getGraph().getEdges()));
    }

}