import org.apache.commons.math3.linear.*;
import org.geotools.graph.structure.*;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by ilias on 20/8/2016.
 */
@SuppressWarnings("unchecked")

public abstract class PhysarumPolycephalum {


    protected Graph graph;
    protected List<Node> sourceNodes;
    protected List<Node> sinkNodes;
    protected final double Io;
    public final double γ;
    protected final int numberOfIterations;
    protected Node sinkNode;
    protected Node sourceNode;

    protected static final Logger logger = Logger.getLogger(PhysarumPolycephalum.class.getName());

    protected int iteration=0;

    public HashMap<Edge, XYSeries> getFlowMap() {
        return flowMap;
    }

    public void setFlowMap(HashMap<Edge, XYSeries> flowMap) {
        this.flowMap = flowMap;
    }

    public HashMap<Edge, XYSeries> getCondMap() {
        return condMap;
    }

    public void setCondMap(HashMap<Edge, XYSeries> condMap) {
        this.condMap = condMap;
    }

    public HashMap<Node, XYSeries> getPressMap() {
        return pressMap;
    }

    public void setPressMap(HashMap<Node, XYSeries> pressMap) {
        this.pressMap = pressMap;
    }

    private HashMap<Edge, XYSeries> flowMap = new HashMap<Edge, XYSeries>();
    private HashMap<Edge, XYSeries> condMap=new HashMap<Edge, XYSeries>();
    private HashMap<Node, XYSeries> pressMap=new HashMap<Node, XYSeries>();

    public PhysarumPolycephalum(Graph graph, List<Node> sourceNodes, List<Node> sinkNodes, double Io, double γ, int numberOfIterations) {
        this.graph = graph;
        this.sourceNodes = sourceNodes;
        this.sinkNodes = sinkNodes;
        this.Io = Io;
        this.γ = γ;
        this.numberOfIterations = numberOfIterations;

    }

    public Graph getGraph() {
        return graph;
    }

    /**
     * Executes the algorithm
     */
    public void execute() {


        logger.info("Starting iterations...");

        for (int i = 0; i <numberOfIterations ; i++) {

            logger.info("Current iteration: "+i);

            sourceNode = chooseRandomNode(sourceNodes);
            do {
                sinkNode = chooseRandomNode(sinkNodes);
            } while (sourceNode.equals(sinkNode));

            List<Node> allNodes = getAllNodes(graph);
            List<Node> allButSinkNodes = getAllButSink(graph);


            RealVector constants = createConstantsVector(allNodes);

            RealMatrix coefficients = createCoefficientsMatrix(allNodes, allButSinkNodes);



           // DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

            DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();
            RealVector solution = solver.solve(constants);

            redefinePressures(allButSinkNodes, solution);

            redefineFlows(graph);

            redefineDiameters(graph);


            iteration++;
        }
    }



    /**
     * @param allNodes All the nodes in the graph
     * @param allButSink All the nodes int the graph, except for the sink node
     * @return returns the coefficients matrix of the linear system
     */


    protected RealMatrix createCoefficientsMatrix(List<Node> allNodes, List<Node> allButSink) {

        double[][] matrix = new double[allNodes.size()][allButSink.size()];

        for (int i = 0; i < allNodes.size(); i++) {
            Node n1 = allNodes.get(i);
            for (int j = 0; j < allButSink.size(); j++) {
                Node n2 = allButSink.get(j);

                if(n1.equals(n2)){
                    matrix[i][j]=calculateSelfCoefficient(n1);
                }else{
                    Edge n1Ton2 =n1.getEdge(n2);
                    if(n1Ton2==null){
                        matrix[i][j]=0;
                    }else{
                        matrix[i][j]=-calculateCoefficient(n1Ton2);
                    }
                }


            }
        }


        return new Array2DRowRealMatrix(matrix);
    }



    /**
     * @param solution Redefines pressures on edges from the solution of the system
     */
    protected void redefinePressures(List<Node> allButSink, RealVector solution) {
        for (int i = 0; i < allButSink.size(); i++) {

            Node n = allButSink.get(i);
            //noinspection Since15
            pressMap.putIfAbsent(n,new XYSeries("("+n.getID()+")"));
            setPressure(n,solution.getEntry(i));
            pressMap.get(n).add(iteration,solution.getEntry(i));
        }


    }

    /**
     *
     * @param graph
     * redefines the flows for the edges of the graph
     */
    @SuppressWarnings("Since15")
    protected void redefineFlows(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        for(Edge e:edges){
            flowMap.putIfAbsent(e,new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double Q =calculateTubeFlow(e);
            setFlow(e,Q);
            flowMap.get(e).add(iteration,Q);
        }

    }


    /**
     *
     * @param graph
     * redefines the flows for the edges of the graph
     */
    protected void redefineDiameters(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        for(Edge e:edges){
            //noinspection Since15
            condMap.putIfAbsent(e, new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double D =calculateTubeDiameter(e);
            setDiameter(e,D);
            condMap.get(e).add(iteration,D);
        }

    }


    /**
     * @param allNodes
     * @return the constants vector of the linear system
     */

    protected RealVector createConstantsVector(List<Node> allNodes) {
        double[] constants = new double[allNodes.size()];
        for (int i = 0; i < constants.length; i++) {

            Node n1 = allNodes.get(i);
            if (n1.equals(sourceNode)) {
                Edge nToSink = n1.getEdge(sinkNode);
                if (nToSink == null) {
                    constants[i] = Io;
                } else {
                    constants[i] = Io + calculateCoefficient(nToSink) * getPressure(sinkNode);
                }

            } else if (n1.equals(sinkNode)) {
                constants[i] = -Io - calculateSelfCoefficient(n1) * getPressure(sinkNode);
            } else {
                Edge nToSink = n1.getEdge(sinkNode);
                if (nToSink == null) {
                    constants[i] = 0;
                } else {
                    constants[i] = calculateCoefficient(nToSink) * getPressure(sinkNode);
                }

            }
        }


        return new ArrayRealVector(constants);
    }

    /**
     * show the pressure diagram
     */
    public void showPressureMap(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );

/*
        Collection<Node> nodes =graph.getNodes();
        for(Node n:nodes){
            dataset.addSeries(pressMap.get(n));
        }
    */
        Iterator it= pressMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }
        JFreeChart flowChart = ChartFactory.createXYLineChart("Time-Pressure", "Time", "H", dataset,PlotOrientation.VERTICAL,true,true,true);
        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");

        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);


    }
    public void showConductivityMap(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );
        Iterator it= condMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }

        JFreeChart flowChart =ChartFactory.createXYLineChart("Time-Conductivity", "Time", "D", dataset, PlotOrientation.VERTICAL,true,false,false);
        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);


    }


    public void showFlowDiagram(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );

        Iterator it= flowMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }
        JFreeChart flowChart =ChartFactory.createXYLineChart("Time-Flow", "Time", "Q", dataset,PlotOrientation.VERTICAL,true,true,true);

        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }




    /**
     * @param nodeList
     * @return A random node from the nodeList
     */

    protected Node chooseRandomNode(List<Node> nodeList) {
        Random random = new Random();


        return nodeList.get(random.nextInt(nodeList.size()));
    }

    protected List<Node> getAllButSink(Graph graph) {
        GraphVisitor allButSinkVisitor = new GraphVisitor() {
            public int visit(Graphable component) {
                if (component.equals(sinkNode)) {
                    return -1;
                } else {
                    return 1;
                }
            }
        };

        return graph.queryNodes(allButSinkVisitor);
    }


    /**
     * @param graph
     * @return Returns a list with all the nodes in the graph
     */

    protected List<Node> getAllNodes(Graph graph) {
        GraphVisitor allNodesVisitor = new GraphVisitor() {
            public int visit(Graphable component) {
                return 1;
            }
        };
        return graph.queryNodes(allNodesVisitor);

    }

    /**
     * @param n
     * @return Returns the self-coefficient
     */

    protected double calculateSelfCoefficient(Node n) {
        double selfC = 0;
        List<Edge> edges = n.getEdges();
        Iterator<Edge> iterator = edges.iterator();
        while (iterator.hasNext()) {
            Edge edge = iterator.next();
            selfC = selfC + calculateCoefficient(edge);

        }


        return selfC;
    }

    /**
     *
     * @param e
     * @return The flow of the edge
     */
    public  abstract double getFlow(Edge e);

    /**
     * @param e
     * @return double
     * <p>
     * Returns the tube length attribute of the edge
     */
    public abstract double getLength(Edge e);

    /**
     * @param e
     * @return Returns the tube diameter attribute of the edge e
     */
    public abstract double getDiameter(Edge e);

    /**
     * @param n
     * @return Returns the pressure attribute of the node n
     */

    public abstract double getPressure(Node n);

    /**
     * @param e
     * @return Calculates  the typical coefficient for the matrix of the system to be solved
     */

    public abstract double calculateCoefficient(Edge e);

    /**
     * @param e
     * @return Recalculates the flow of the tube according to the formula used in the problem
     */


    public abstract double calculateTubeFlow(Edge e);

    /**
     * @param e
     * @return Recalculates the diameter of the tube according to the sigmoidal response
     * of the diameter
     */

    public abstract double calculateTubeDiameter(Edge e);

    /**
     * @param n sets the pressure for the node
     */
    public abstract void setPressure(Node n, double p);

    /**
     * @param e sets the flow for the edge
     */
    public abstract void setFlow(Edge e, double Q);

    /**
     * @param e sets the diameter of the edge
     */
    public abstract void setDiameter(Edge e, double D);
}
