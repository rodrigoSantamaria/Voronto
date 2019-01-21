package es.usal.voronto.view;

import javax.swing.JFrame;
import javax.swing.JTextField;

import es.usal.voronto.control.Voronto;

import java.awt.Color;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.JLabel;

public class SearchFrame extends JFrame implements KeyListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = -760578346385933961L;
	public JTextField textField;
	private VoronoiVisualization vv;
	public SearchFrame(Voronto v) {
		this.setBounds(0, 0, 220, 50);
		getContentPane().setBackground(Color.WHITE);
		setBackground(Color.WHITE);
		setAlwaysOnTop(true);
		setResizable(false);
		getContentPane().setLayout(null);
		this.setUndecorated(true);
		
		//this.setContentPane(v);
		this.vv=v.gv;
		
		textField = new JTextField();
		textField.setToolTipText("Write here any gene or term name");
		textField.setBounds(6, 21, 208, 26);
		getContentPane().add(textField);
		textField.setColumns(10);
		textField.addKeyListener(this);
		
		JLabel lblSearchFor = new JLabel("Search for:");
		lblSearchFor.setForeground(Color.GRAY);
		lblSearchFor.setBackground(Color.WHITE);
		lblSearchFor.setBounds(73, 6, 87, 16);
		getContentPane().add(lblSearchFor);
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
		switch(e.getKeyCode())
			{
			case KeyEvent.VK_ENTER:
				vv.search(textField.getText());
				break;
			case KeyEvent.VK_ESCAPE:
				vv.search(null);
				break;
			}
	}
}
