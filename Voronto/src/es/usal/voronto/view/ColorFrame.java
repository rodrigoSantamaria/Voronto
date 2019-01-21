package es.usal.voronto.view;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JTextField;

import es.usal.voronto.control.Voronto;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.JLabel;
import javax.swing.JPanel;

public class ColorFrame extends JFrame implements KeyListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = -760578346385933961L;
	
	private VoronoiVisualization vv;
	public ColorChooser[] choosers;
	/**
	 * Returns a <code>JPanel</code> with the interface for the configuration of diagram colors
	 * 
	 * @param paleta Array with the diagram colors
	 * @param textoLabel Array with labels that define each color
	 * @param muestra <code>JTextField</code> where the colors to be shown are to be drawn
	 * @return <code>JPanel</code> with the configuration interface of diagram colors
	 */
	public ColorFrame(Voronto v, Color[] paleta, String[] textoLabel, JTextField[] muestra) {
		this.setBounds(0, 0, 280, 180);
		getContentPane().setBackground(Color.WHITE);
		setBackground(Color.WHITE);
		//setAlwaysOnTop(true);
		setResizable(false);
		getContentPane().setLayout(null);
		this.setUndecorated(true);
		
		//this.setContentPane(v);
		this.vv=v.gv;
		
		
		
		//Configuramos el panel de configuracion de la paleta de colores
		JPanel panelColor = new JPanel();
		//JFrame panelColor=new JFrame();
		panelColor.setLayout(new GridBagLayout());
		
		JLabel[] etiqueta = new JLabel[paleta.length];
		JButton[] boton = new JButton[paleta.length];
		choosers=new ColorChooser[paleta.length];
		
		GridBagConstraints constraints = new GridBagConstraints();
		for(int i = 0; i < paleta.length; i++){
			boton[i] = new JButton();
			etiqueta[i] = new JLabel(textoLabel[i]);
			muestra[i] = new JTextField(10);
			muestra[i].setEditable(false);
			muestra[i].setBackground(paleta[i]);
			boton[i].setText("Change");
			choosers[i]=new ColorChooser(boton[i],muestra[i], paleta[i]);
			boton[i].addActionListener(choosers[i]);
			
			constraints.gridx = 0;
			constraints.gridy = i;
			constraints.anchor = GridBagConstraints.WEST;
			panelColor.add(etiqueta[i],constraints);
			
			constraints.anchor = GridBagConstraints.CENTER;
			constraints.gridx = 1;
			constraints.fill = GridBagConstraints.HORIZONTAL;
			constraints.weightx = 1.0;
			panelColor.add(muestra[i],constraints);
			
			constraints.gridx = 2;
			constraints.weightx = 0.0;
			panelColor.add(boton[i],constraints);
		}
		
		
		
		panelColor.setBackground(Color.WHITE);
		panelColor.setBounds(6, 6, 260, 116);
		getContentPane().add(panelColor);
		
		
		JButton cancelButton = new JButton("Done");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				setVisible(false);
				Color[] colors=new Color[3];
				int i=0;
				for(ColorChooser c:choosers)
						colors[i++]=c.color;
				vv.changeColors(colors);
				}
			});
		cancelButton.setBounds(100, 130, 80, 30);
		cancelButton.setActionCommand("Done");
		cancelButton.setVisible(true);
		getContentPane().add(cancelButton);
		
	}
	
	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void keyReleased(KeyEvent e) {
		switch(e.getKeyChar())
			{
			case 'c':
			case 'q':
			case 'e':
				this.setVisible(false);
				Color[] colors=new Color[3];
				int i=0;
				for(ColorChooser c:choosers)
					colors[i++]=c.color;
				vv.changeColors(colors);
				
			}
	/*	switch(e.getKeyCode())
			{
			case KeyEvent.VK_ENTER:
				vv.search(textField.getText());
				break;
			case KeyEvent.VK_ESCAPE:
				vv.search(null);
				break;
			}*/
	}
}
