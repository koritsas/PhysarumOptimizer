package edu.koritsas.slimemold;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;
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
import java.util.stream.Collectors;

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
    private HashMap<Edge,Double> fluxMap ;
    private HashMap<Edge,Double> conductivityMap;
    private HashMap<Node,Double> pressureMap;

    protected static final Logger logger = Logger.getLogger(PhysarumPolycephalum.class.getName());

    protected int iteration=0;

    private HashMap<Edge, XYSeries> flowSeriesMap = new HashMap<Edge, XYSeries>();
    private HashMap<Edge, XYSeries> conductivitySeriesMap =new HashMap<Edge, XYSeries>();
    private HashMap<Node, XYSeries> pressureSeriesMap =new HashMap<Node, XYSeries>();

    public HashMap<Edge, XYSeries> getFlowSeriesMap() {
        return flowSeriesMap;
    }

    public void setFlowSeriesMap(HashMap<Edge, XYSeries> flowSeriesMap) {
        this.flowSeriesMap = flowSeriesMap;
    }

    public HashMap<Edge, XYSeries> getConductivitySeriesMap() {
        return conductivitySeriesMap;
    }

    public void setConductivitySeriesMap(HashMap<Edge, XYSeries> conductivitySeriesMap) {
        this.conductivitySeriesMap = conductivitySeriesMap;
    }

    public HashMap<Node, XYSeries> getPressureSeriesMap() {
        return pressureSeriesMap;
    }

    public void setPressureSeriesMap(HashMap<Node, XYSeries> pressureSeriesMap) {
        this.pressureSeriesMap = pressureSeriesMap;
    }


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

        initializeMaps(graph);

        logger.info("Starting iterations...");

        for (int i = 0; i <numberOfIterations ; i++) {

            logger.info("Current iteration: "+i);
            List<Node> allNodes = getAllNodes(graph);

            chooseRandomNodes(allNodes);




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

        eliminateEdges(graph);
    }

    private void chooseRandomNodes(List<Node> allNodes) {
       sourceNode = chooseRandomNode(sourceNodes);
        do {
            sinkNode = chooseRandomNode(sinkNodes);
        } while (sourceNode.equals(sinkNode));





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

    private void initializeMaps(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        fluxMap=new HashMap<Edge, Double>(graph.getEdges().size());
        conductivityMap=new HashMap<Edge,Double>(graph.getEdges().size());
        for (Edge e:edges){
            fluxMap.putIfAbsent(e,0.0);
            conductivityMap.putIfAbsent(e,1.0);
        }

        Collection<Node> nodes =graph.getNodes();
        pressureMap=new HashMap<Node,Double>(nodes.size());
        for (Node n:nodes){
            pressureMap.putIfAbsent(n,0.0);

        }




    }



    /**
     * @param solution Redefines pressures on edges from the solution of the system
     */
    protected void redefinePressures(List<Node> allButSink, RealVector solution) {
        for (int i = 0; i < allButSink.size(); i++) {

            Node n = allButSink.get(i);
            //noinspection Since15
            pressureSeriesMap.putIfAbsent(n,new XYSeries("("+n.getID()+")"));
           pressureMap.put(n,solution.getEntry(i));
            pressureSeriesMap.get(n).add(iteration,solution.getEntry(i));
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
            flowSeriesMap.putIfAbsent(e,new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double Q =calculateTubeFlux(e);
            fluxMap.put(e,Q);
            flowSeriesMap.get(e).add(iteration,Q);
        }

    }

    public double getSolutionCost() {

        Collection<Edge> edges = graph.getEdges();
        double cost = edges.stream().collect(Collectors.summingDouble(value -> getEdgeCost(value)));
        return  cost;
    }
    private void eliminateEdges(Graph graph){

        graph.getEdges().removeIf(o -> FastMath.abs(conductivityMap.get(o))<0.0001);


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
            conductivitySeriesMap.putIfAbsent(e, new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double D =calculateTubeConductivity(e);
            conductivityMap.put(e,D);
            conductivitySeriesMap.get(e).add(iteration,D);
        }

    }

    private double calculateTubeFlux(Edge e){

        double p1 =pressureMap.get(e.getNodeA());
        double p2 =pressureMap.get(e.getNodeB());

        double D =conductivityMap.get(e);
        double w = getEdgeCost(e);

        double Q=(D/w)*(p1-p2);

        return Q;
    }

    private double calculateTubeConductivity(Edge e){
        double D = conductivityMap.get(e);

        double Q = fluxMap.get(e);
        double L=getEdgeCost(e);

        double p1 =pressureMap.get(e.getNodeA());
        double p2 =pressureMap.get(e.getNodeB());

        double ps=pressureMap.get(sourceNode);
        double pe =pressureMap.get(sinkNode);

       // double newD =(0.5)*((Q*(p1-p2))/(L*(ps-pe))+D);

        double fQ = Math.pow(Math.abs(Q), γ) / (1 + Math.pow(Math.abs(Q), γ));
        //double fQ=Math.abs(Q);

        double newD = fQ - 0.1 * D;



        return newD;


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
                    constants[i] = Io + calculateCoefficient(nToSink) * pressureMap.get(sinkNode);
                }

            } else if (n1.equals(sinkNode)) {
                constants[i] = -Io - calculateSelfCoefficient(n1) *  pressureMap.get(sinkNode);
            } else {
                Edge nToSink = n1.getEdge(sinkNode);
                if (nToSink == null) {
                    constants[i] = 0;
                } else {
                    constants[i] = calculateCoefficient(nToSink) *  pressureMap.get(sinkNode);
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
            dataset.addSeries(pressureSeriesMap.get(n));
        }
    */
        Iterator it= pressureSeriesMap.entrySet().iterator();
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
        Iterator it= conductivitySeriesMap.entrySet().iterator();
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

        Iterator it= flowSeriesMap.entrySet().iterator();
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

    private double calculateCoefficient(Edge e){
        double D =conductivityMap.get(e);
        double w = getEdgeCost(e);

        return D/w;
    }

    /**
     *
     * @param e
     * @return the weight of the edge e
     */
    public abstract double getEdgeCost(Edge e);

    /**
     *
     * @param e
     * @return the constraint value of the edge e
     */
    public abstract double getEdgeConstraint(Edge e);


    public abstract double getEuclideanDistance(Node n1,Node n2);


}
